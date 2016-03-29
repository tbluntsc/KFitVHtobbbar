import os
import ROOT
import numpy as np
from numpy.linalg import inv

#Note the _small in the name of the root file: to try out new features
#loading all 25 files takes too long, so a smaller version of the root
#file will be created, only reading one of the root file in "address"
out = ROOT.TFile("ChiSquareFits.root", "UPDATE")
address = "dcap://t3se01.psi.ch:22125////pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-higgs/run2/VHBBHeppyV20/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/VHBB_HEPPY_V20_ZH_HToBB_ZToLL_M125_13TeV_powheg_Py8__fall15MAv2-pu25ns15v1_76r2as_v12-v1/160209_172236/0000/"

#Create an empty Tree 
tree = ROOT.TTree("ChiSquareFits", "Tree containing results from kinematic fit")

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

print "Total number of entries: ", chain.GetEntries()
counter = 0.0
#To speed up runtime 3e+3 instead of 1e+11
for iev in range(int( min(3e+3, chain.GetEntries()))):
    chain.GetEntry(iev)
    ev = chain
    if iev%500 == 0:
        print "Processing event ", iev+1
    
    #Discard all entries w Vtype not equal to either 0 or 1. (1: V -> e+e- / e-v_e, 0: V -> mumu / mu v_mu )
    if (ev.Vtype != 0) & (ev.Vtype != 1):
        print "VType neither 0 or 1 in event ", str(iev), " Vtype = ", str(ev.Vtype) 
        continue

    #Discard all entries w Vectorboson Momentum smaller than 50 GeV
    if ev.V_pt < 50:
        print "V_pt < 50 in event ", str(iev), " V_pt = ", str(ev.V_pt)
        continue

    #Discard all entries w Higgs-tagged Jets that have PT < 20 or |eta| > 2.4 or hadronFlavour != 5
    higgs_bools = []
    for H_jets in xrange(len(ev.hJidx)):
        if (ev.Jet_pt[ev.hJidx[H_jets]] < 20) or (ev.Jet_eta[ev.hJidx[H_jets]] > 2.4) or (ev.Jet_eta[ev.hJidx[H_jets]] < -2.4) or (ev.Jet_hadronFlavour[ev.hJidx[H_jets]] != 5):
            higgs_bools.append(True)
    if any(higgs_bools):
        print "Higgs jet w pt < 20, |eta| > 2.4 or Flavour != 5 in event ", str(iev)
        continue

    # For each event we want to build three matrices: 
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

    # The second matrix is the covariance matrix. In our case this matrix is simply a diagonal matrix w. the variances squared as entries
    # V = diag(sigma^2(E_1),..,sigma^2(E_nJet))
    # To get these sigmas we use the root file Sigmas_Fit.root where there are histograms w. a plot of sigma as a function of the Energy for different
    # eta regions and hadronFlavours. For each event we thus extract Eta, hadronFlavour & the energies of the jets and then fill a diagonal matrix
    # w the fit from the root file.
    # eta0: eta in [0.0,1.0]
    # eta1: eta in [1.0,1.5]
    # eta2: eta in [1.5,2.0]
    # eta3: eta in [2.0,2.5]
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
        print "Jet w Flavour not in [0,5] in event " , str(iev)
        continue
        
#    Sigma_root = ROOT.TFile("Sigmas_Fit.root")
    regions = np.zeros(ev.nJet)
    for jets in xrange(len(jet_etas)):
        if (jet_etas[jets] > 0.0 and jet_etas[jets] < 1.0):
            regions[jets] = 0
        elif (jet_etas[jets] > 1.0 and jet_etas[jets] < 1.5):
            regions[jets] = 1
        elif (jet_etas[jets] > 1.5 and jet_etas[jets] < 2.0):
            regions[jets] = 2
        elif (jet_etas[jets] > 2.0 and jet_etas[jets] < 2.5):
            regions[jets] = 3
        else: 
            regions[jets] = 99
    if any(x == 99 for x in regions):
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

    # The last matrix we need to build is the matrix L coming from the constraints, defined by L dot Theta = R 
    # R = [(- cos (V_phi) * V_pt),(-sin(V_phi)* V_pt)]
    # L is a 2 x nJet matrix, w. the first row containing the cosines of the Jet_phis, the second row containing the sines of the Jet_phis

    R = np.zeros(2)
    R[0] = -np.cos(ev.V_phi)*ev.V_pt
    R[1] = -np.sin(ev.V_phi)*ev.V_pt

    first_row = np.zeros(ev.nJet)
    second_row = np.zeros(ev.nJet)

    for i in xrange(ev.nJet):
        first_row[i] = np.cos(ev.Jet_phi[i])
        second_row[i] = np.sin(ev.Jet_phi[i])

    L = np.matrix([first_row, second_row])

    # The only thing left to do is calculate the various matrix products.

    A_tr = np.transpose(A)
    V_inv = inv(V)
    L_tr = np.transpose(L)

    C = np.dot(A_tr, np.dot(V_inv, A))
    C_inv = inv(C)
    
    LC1LT1 = np.dot(L,np.dot(C_inv, L_tr))
    LC1LT1_inv = inv(LC1LT1)
    F = C_inv - np.dot(C_inv,np.dot(L_tr,np.dot(LC1LT1_inv,np.dot(L, C_inv))))
    G = np.dot(LC1LT1_inv,np.dot(L,C_inv))
    H = -1.0*LC1LT1_inv
    
    Theta = np.dot(np.dot(F, np.dot(A_tr,V_inv)),jet_pts) + np.dot(np.transpose(G),R)
    print "Estimates :",Theta
    print "Measured values :", jet_pts

    counter += 1

print "We found ",counter, " fitting events"    
a = raw_input("We did it boys")
out.cd()
tree.Print()
tree.Scan("GenPt")
tree.Scan("nJet")
tree.Scan("hJidx")
tree.Write("", ROOT.TObject.kOverwrite)
out.Close()
