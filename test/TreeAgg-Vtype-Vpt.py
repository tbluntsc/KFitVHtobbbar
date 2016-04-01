import os
import ROOT
import numpy as np
from numpy.linalg import inv

print_discriminating_reasons = False

#Create a new ROOT File where the Chi square fits will be saved in the end & load the address of the data
out = ROOT.TFile("ChiSquareFits.root", "UPDATE")
address = "dcap://t3se01.psi.ch:22125////pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-higgs/run2/VHBBHeppyV20/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/VHBB_HEPPY_V20_ZH_HToBB_ZToLL_M125_13TeV_powheg_Py8__fall15MAv2-pu25ns15v1_76r2as_v12-v1/160209_172236/0000/"

#Create an empty Tree 
tree = ROOT.TTree("ChiSquareFits", "Tree containing results from kinematic fit")

#Initialize branches of tree
maxJet=100
nJet = np.zeros(1, dtype = int)
estimates = np.zeros(maxJet, dtype = np.float64())
measurements = np.zeros(maxJet, dtype = np.float64())
etas = np.zeros(maxJet, dtype = np.float64())
higgs_tag = np.zeros(maxJet, dtype = int)
est_mass = np.zeros(1, dtype = np.float64())
mea_mass = np.zeros(1, dtype = np.float64())

tree.Branch('nJet', nJet, 'nJet/I')
tree.Branch('estimates', estimates, "estimates[nJet]/D")
tree.Branch('measurements', measurements, "measurements[nJet]/D")
tree.Branch('etas', etas, "etas[nJet]/D")
tree.Branch('higgs_tag', higgs_tag, "higgs_tag[nJet]/I")
tree.Branch('est_mass', est_mass, "est_mass/D")
tree.Branch('mea_mass', mea_mass, "mea_mass/D")

#Create array w. strings corresponding to tree names (i.e. tree_1.root, tree_2.root etc.)
no_of_files = 1
file_names = []
for i in range(1, no_of_files+1):
    file_names.append("tree_"+str(i)+".root")

print "Added ", len(file_names), "files"

#Loop over trees to chain them together
chain = ROOT.TChain("tree")
CountWeighted = 0

for file_name in file_names:
    f = ROOT.TFile.Open(address + "/" + file_name)
    if f==None or f.IsZombie():
        print "File/Address does not contain a tree"
        continue
    f.cd()
    count = f.Get("CountWeighted")
    CountWeighted += count.GetBinContent(1)
    f.Close()
    chain.AddFile(address + "/" + file_name)

print "Total processed events: ", CountWeighted

#Initialize all the branches needed for the computation
chain.SetBranchStatus("*", False)
chain.SetBranchStatus("Jet_hadronFlavour", True)
chain.SetBranchStatus("Jet_mcPt", True)
chain.SetBranchStatus("Jet_pt", True)
chain.SetBranchStatus("Jet_eta", True)
chain.SetBranchStatus("Jet_phi", True)
chain.SetBranchStatus("Jet_mass", True)
chain.SetBranchStatus("Jet_mcEta", True)
chain.SetBranchStatus("Jet_mcPhi", True)
chain.SetBranchStatus("Jet_pt_reg", True)
chain.SetBranchStatus("nJet", True)
chain.SetBranchStatus("Vtype", True)
chain.SetBranchStatus("V_pt", True)
chain.SetBranchStatus("hJidx", True)
chain.SetBranchStatus("V_phi", True)

print "Total number of entries: ", chain.GetEntries()

#Initialize counter that conts the number of events for which fitted jet energies were calculated and the event was not discarded.
counter = 0.0

#To speed up runtime 3e+3 instead of 1e+11, i.e. only load the first 3000 events.
for iev in range(int( min(2e+4, chain.GetEntries()))):
    chain.GetEntry(iev)
    ev = chain
    if iev%500 == 0:
        print "Processing event ", iev+1

    #Discard all entries w Vtype not equal to either 0 or 1. (1: V -> e+e- / e-v_e, 0: V -> mumu / mu v_mu )
    if (ev.Vtype != 0) & (ev.Vtype != 1):
        if print_discriminating_reasons:
            print "VType neither 0 or 1 in event ", str(iev), " Vtype = ", str(ev.Vtype) 
        continue

    #Discard all entries w Vectorboson Momentum smaller than 50 GeV
    if ev.V_pt < 50:
        if print_discriminating_reasons:
            print "V_pt < 50 in event ", str(iev), " V_pt = ", str(ev.V_pt)
        continue

    #Discard all entries w Higgs-tagged Jets that have PT < 20 or |eta| > 2.4 or hadronFlavour != 5
    higgs_bools = []
    for H_jets in xrange(len(ev.hJidx)):
        if (ev.Jet_pt[ev.hJidx[H_jets]] < 20) or (ev.Jet_eta[ev.hJidx[H_jets]] > 2.4) or (ev.Jet_eta[ev.hJidx[H_jets]] < -2.4) or (ev.Jet_hadronFlavour[ev.hJidx[H_jets]] != 5):
            higgs_bools.append(True)
    if any(higgs_bools):
        if print_discriminating_reasons:
            print "Higgs jet w pt < 20, |eta| > 2.4 or Flavour != 5 in event ", str(iev)
        continue

    #Discard events where only one jet was produced (Is this unphysical? Or what's the problem? LC1LT seems to always be singular in that case)
    if ev.nJet == 1:
        if print_discriminating_reasons:
            print "Only one jet in event ", iev
        continue

#    if ev.nJet == 2:
#        if print_discriminating_reasons:
#            print "Only two jets in event ", iev
#        continue
    # For each event we want to build three matrices, A (model matrix), V (covariance matrix) & L (constraint matrix)
    #
    # The model matrix A. This matrix tells us in what way the measured values and the estimated values are related. It is defined by
    # Y = A*Theta + epsilon where
    # Y: measured values
    # Theta: values we want to estimate
    # epsilon: offset
    #
    # in this case we assume a linear detector, i.e. Y = a*Theta + b with a being a scalar and b being a vector filled w the scalar b.
    # The dimension of this matrix is n x n where n is the number of jets we want to estimate (i.e. nJet)
    
    a = 1.0
    b = 0.0
    diagonal = np.zeros(ev.nJet)
    for jets in xrange((ev.nJet)):
        diagonal[jets] = a
    A = np.diag(diagonal)
    if np.linalg.matrix_rank(A) != A.shape[0]:
        print "Matrix A is not of full rank in event ", iev
        continue

    # The second matrix is the covariance matrix. In our case this matrix is simply a diagonal matrix w. the variances squared as entries
    # V = diag(sigma^2(E_1),..,sigma^2(E_nJet))
    # To get these sigmas we use the root file Sigmas_Fit.root where there are histograms w. a plot of sigma as a function of the Energy for different
    # eta regions and hadronFlavours. For each event we thus extract Eta, hadronFlavour & the energies of the jets and then fill a diagonal matrix
    # w the fit from the root file.
    # eta0: abs(eta) in [0.0,1.0]
    # eta1: abs(eta) in [1.0,1.5]
    # eta2: abs(eta) in [1.5,2.0]
    # eta3: abs(eta) in [2.0,2.5]

    jet_pts = np.zeros(ev.nJet)
    jet_sigmas = np.zeros(ev.nJet)
    jet_etas = np.zeros(ev.nJet)
    jet_flavours = np.zeros(ev.nJet)
    jet_region = np.zeros(ev.nJet)

    for jets in xrange(ev.nJet):
        jet_etas[jets] = ev.Jet_eta[jets]
        jet_flavours[jets] = ev.Jet_hadronFlavour[jets]
        jet_pts[jets] = ev.Jet_pt[jets]

    #If any of the flavours of the jets are not 0 or 5, go to the next event. This is because sigmas were only fitted for Flavour = 0 or Flavour = 5
    if any(x not in [0,5] for x in jet_flavours):
        if print_discriminating_reasons:
            print "Jet w Flavour not in [0,5] in event " , str(iev)
        continue
        
    regions = np.zeros(ev.nJet)
    for jets in xrange(len(jet_etas)):
        if (np.absolute(jet_etas[jets]) > 0.0 and np.absolute(jet_etas[jets]) < 1.0):
            regions[jets] = 0
        elif (np.absolute(jet_etas[jets]) > 1.0 and np.absolute(jet_etas[jets]) < 1.5):
            regions[jets] = 1
        elif (np.absolute(jet_etas[jets]) > 1.5 and np.absolute(jet_etas[jets]) < 2.0):
            regions[jets] = 2
        elif (np.absolute(jet_etas[jets]) > 2.0 and np.absolute(jet_etas[jets]) < 2.5):
            regions[jets] = 3
        else: 
            regions[jets] = 99
    if any(x == 99 for x in regions):
        if print_discriminating_reasons:
            print "Jet w Eta not in any of the 4 regions in event ", str(iev), regions, jet_etas
        continue


    histo_strings = []
    for jets in xrange(len(regions)):

        string = []
        string.append("Sigmas_")
        string.append(str(int(regions[jets])))
        string.append("_")
        string.append(str(int(jet_flavours[jets])))
        string = "".join(string)
        histo_strings.append(string)

    for jets in xrange(len(jet_sigmas)):

        RootFile = ROOT.TFile("Sigmas_Fits.root")
        RootFile.cd()

        current_histo = ROOT.gDirectory.Get(histo_strings[jets])
        myfunc = current_histo.GetFunction("sigma_func")

        jet_sigmas[jets] = myfunc.Eval(ev.Jet_pt[jets])

    diagonal = np.zeros(ev.nJet)
    for jets in xrange((ev.nJet)):
        diagonal[jets] = jet_sigmas[jets]**2
    V = np.diag(diagonal)
    if np.linalg.matrix_rank(V) !=V.shape[0]:
        print "Matrix V is singular in event ", iev
        continue

    # The last matrix we need to build is the matrix L coming from the constraints, defined by L dot Theta = R (Theta being the vector containing the measurements)
    # R = [(- cos (V_phi) * V_pt),(-sin(V_phi)* V_pt)]
    # L is a 2 x nJet matrix, w. the first row containing the sines of the Jet_phis, the second row containing the cosines of the Jet_phis

    Lorentzvectors = []
    for i in xrange(ev.nJet):
        v = ROOT.TLorentzVector()
        v.SetPtEtaPhiM(ev.Jet_pt[i], ev.Jet_eta[i], ev.Jet_phi[i], ev.Jet_mass[i])
        Lorentzvectors.append(v)

    lepton_vector = ROOT.TLorentzVector()
    lepton_vector.SetPtEtaPhiM(ev.V_pt, ev.V_eta, ev.V_phi, ev.V_mass)
    Lorentzvectors.append(lepton_vector)

    R = np.zeros(2)
    R[0] = -lepton_vector.Px()
    R[1] = -lepton_vector.Py()

    jet_phis = np.zeros(ev.nJet)
    for i in xrange(ev.nJet):
        jet_phis[i] = ev.Jet_phi[i]

    L = np.matrix([np.cos(jet_phis), np.sin(jet_phis)])
 
        
   # The only thing left to do is calculate the various matrix products.

    A_tr = np.transpose(A)
    V_inv = inv(V)
    L_tr = np.transpose(L)

    C = np.dot(A_tr, np.dot(V_inv, A))
    if np.linalg.matrix_rank(C) !=C.shape[0]:
        print "C matrix not of full rank in event ", iev
        continue
    C_inv = inv(C)
    
    LC1LT1 = np.dot(L,np.dot(C_inv, L_tr))
    if np.linalg.matrix_rank(LC1LT1) != LC1LT1.shape[0]:
        print "LC1LT1 matrix not of full rank in event ", iev
        continue

    LC1LT1_inv = inv(LC1LT1)
    F = C_inv - np.dot(C_inv,np.dot(L_tr,np.dot(LC1LT1_inv,np.dot(L, C_inv))))
    G = np.dot(LC1LT1_inv,np.dot(L,C_inv))
    H = -1.0*LC1LT1_inv

    #Theta are the estimated values for the jet pts
    Theta = np.dot(np.dot(F, np.dot(A_tr,V_inv)),jet_pts) + np.dot(np.transpose(G),R)

    Lorentzvectors = []
    for i in xrange(ev.nJet):
        v = ROOT.TLorentzVector()
        v.SetPtEtaPhiM(Theta[0,i],ev.Jet_eta[i], ev.Jet_phi[i], ev.Jet_mass[i])
        Lorentzvectors.append(v)
    Lorentzvectors.append(lepton_vector)
    px = 0.0
    py = 0.0

    for i in xrange(len(Lorentzvectors)):
        px += Lorentzvectors[i].Px()
        py += Lorentzvectors[i].Py()
    print "--------------------------------------------"
    print " "
    print "Fitted jets in event ", iev
    print "Is pT conservation satisfied?"
    print " "
    print "px = ", px
    print "py = ", py
    print " "

    estimat = []
    measure = []
    for i in xrange(len(ev.hJidx)):
        estimat.append(Theta[0,ev.hJidx[i]])
        measure.append(jet_pts[ev.hJidx[i]])

    print "Estimates for Higgs PTs:", estimat
    print "Measured values for Higgs PTs:", measure

    if any(x < 0.0 for x in estimat):
        print "ESTIMATED VALUE NEGATIVE!"

    Higgs_lorentz = []
    for i in xrange(len(ev.hJidx)):
        v = ROOT.TLorentzVector()
        v.SetPtEtaPhiM(estimat[i], jet_etas[ev.hJidx[i]], ev.Jet_phi[ev.hJidx[i]], ev.Jet_mass[ev.hJidx[i]])
        Higgs_lorentz.append(v)
    
    higgs_vector = ROOT.TLorentzVector()
    for i in xrange(len(Higgs_lorentz)):
        higgs_vector += Higgs_lorentz[i]

    Higgs_measured = []
    for i in xrange(len(ev.hJidx)):
        v = ROOT.TLorentzVector()
        v.SetPtEtaPhiM(jet_pts[ev.hJidx[i]], jet_etas[ev.hJidx[i]], ev.Jet_phi[ev.hJidx[i]], ev.Jet_mass[ev.hJidx[i]])
        Higgs_measured.append(v)

    higgs_vector_m = ROOT.TLorentzVector()
    for i in xrange(len(Higgs_measured)):
        higgs_vector_m += Higgs_measured[i]

    print "Estimated Higgs mass: ", higgs_vector.M()
    print "Measured Higgs mass: ", higgs_vector_m.M()
    counter += 1

    for i in xrange(ev.nJet):
        if Theta[0,i] < 0:
            print "theta : ",Theta[0,:]
            print "measurements: ",jet_pts
            print "etas: ", jet_etas
            phis = []
            for j in xrange(ev.nJet):
                phis.append(ev.Jet_phi[j])
            print "phis: ", phis
            print "V_phi: ", ev.V_phi
            print "V_pt: ", ev.V_pt
            #ha = raw_input( "Found negative Jet PT" )

    print Theta.shape
    print ev.nJet
    print Theta[0,ev.nJet-1]
    nJet[0] = ev.nJet
    est_mass[0] = higgs_vector.M()
    mea_mass[0] = higgs_vector_m.M()
    
    for jet in xrange(ev.nJet):
        estimates[jet] = Theta[0,jet]
        measurements[jet] = ev.Jet_pt[jet]
        etas[jet] = ev.Jet_eta[jet]        
        if jet in ev.hJidx:
            higgs_tag[jet] = 1
        else:
            higgs_tag[jet] = 0
    print estimates
    print measurements
    print etas
    print higgs_tag
    tree.Fill()

print "We found ",counter, " fitting events"    

out.cd()
tree.Write("", ROOT.TObject.kOverwrite)
out.Close()
