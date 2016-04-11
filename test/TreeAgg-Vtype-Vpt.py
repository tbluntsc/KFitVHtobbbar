import ROOT
import numpy as np
from numpy.linalg import inv

print_discriminating_reasons = False

# Opening RootFile containing the Sigmas Fits
RootFile = ROOT.TFile("Sigmas_Fits.root")

#Later in the code we will need to create three matrices to solve the Lagrangian system, A (model matrix), V (covariance matrix) & L (constraints matrix)
#Since these matrices will have to be built in each iteration of the loop over events a function is defined for each matrix which can be called inside the loop

def create_cut_permutations(Theta):
    result = []
    for i in xrange(len(Theta)):
        idx = []
        for j in xrange(len(Theta)):
            if j != i:
                idx.append(j)
        result.append(idx)
    return result

def ChiSquare(V, estimated_pt, original_pt):
    V_diag = np.diag(V).tolist()
    res = 0.0
    for idx in xrange(len(estimated_pt)):
        res += ((estimated_pt[idx] - original_pt[idx])**2)/(V_diag[idx])
    return res

def find_minimal_array(Theta, original_pt, V):
    idx_list = create_cut_permutations(np.arange(len(Theta)))

    first_V = V[np.ix_(idx_list[0], idx_list[0])]
    min_chi = ChiSquare(first_V, [Theta[i] for i in idx_list[0]], [original_pt[i] for i in idx_list[0]])
    best_index = idx_list[0]

    for indices in idx_list:
        current_chi = ChiSquare(V[np.ix_(indices, indices)], [Theta[i] for i in indices], [original_pt[i] for i in indices])
        if current_chi < min_chi:
            best_index = indices
            min_chi = current_chi

    return best_index

def A_matrix(a,b, nJet):
    diagonal = np.zeros(nJet)
    for jets in xrange(nJet):
        diagonal[jets] = a
    A = np.diag(diagonal)
    return A

def V_matrix(regions, nJet, jet_pts, jet_etas, jet_flavours, SigmasFile):

    histo_strings = []
    for jet in xrange(len(regions)):
        string = "Sigmas_" + str(int(regions[jet])) + "_" + str(int(jet_flavours[jet]))
        histo_strings.append(string)

    jet_sigmas = np.zeros(nJet)
    for jet in xrange(len(jet_sigmas)):
        SigmasFile.cd()
        current_histo = ROOT.gDirectory.Get(histo_strings[jet])
        myfunc = current_histo.GetFunction("sigma_func")
        
        jet_sigmas[jet] = myfunc.Eval(jet_pts[jet])

    diagonal = np.zeros(nJet)
    for jet in xrange(nJet):
        diagonal[jet] = jet_sigmas[jet]**2
    
    V = np.diag(diagonal)

    return V

def L_matrix_and_R_vector(nJet, jet_pts, V_pt, V_eta, V_phi, V_mass, jet_phis, jet_mass, jet_etas):
    
    Lorentzvectors = []
    for jet in xrange(nJet):
        v = ROOT.TLorentzVector()
        v.SetPtEtaPhiM(jet_pts[jet], jet_etas[jet], jet_phis[jet], jet_mass[jet])
        Lorentzvectors.append(v)

    lepton_vector = ROOT.TLorentzVector()
    lepton_vector.SetPtEtaPhiM(V_pt, V_eta, V_phi, V_mass)

    R = np.zeros(2)
    R[0] = -lepton_vector.Px()
    R[1] = -lepton_vector.Py()

    L = np.matrix([np.cos(jet_phis), np.sin(jet_phis)])

    return L, R

def LagrangianSolver(A, L, V,  R, jet_pts):
    A_tr = np.transpose(A)
    V_inv = inv(V)
    L_tr = np.transpose(L)

    C = np.dot(A_tr, np.dot(V_inv,A))
    C_inv = inv(C)

    LC1LT1 = np.dot(L, np.dot(C_inv, L_tr))
    LC1LT1_inv = inv(LC1LT1)
    F = C_inv - np.dot(C_inv, np.dot(L_tr, np.dot(LC1LT1_inv, np.dot(L, C_inv))))
    G = np.dot(LC1LT1_inv, np.dot(L,C_inv))
    H = -1.0*LC1LT1_inv

    return np.dot(np.dot(F,np.dot(A_tr,V_inv)), jet_pts) + np.dot(np.transpose(G),R)

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
no_of_files = 10
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
chain.SetBranchStatus("V_eta", True)
chain.SetBranchStatus("V_mass", True)

print "Total number of entries: ", chain.GetEntries()

no_fitted_events = 0.0

#To speed up runtime 3e+3 instead of 1e+11, i.e. only load the first 3000 events.
for iev in range(int( min(1e+5, chain.GetEntries()))):
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

    #Discard all events w. only 2 Jets
#    if ev.nJet==2:
#        if print_discriminating_reasons:
#            print "Only two jets in event", iev
#        continue
##### Building matrices & solving Lagrangian system #####

    a = 1.0
    b = 0.0

    A = A_matrix(a,b, ev.nJet)

    if np.linalg.matrix_rank(A) != A.shape[0]:
        print "Matrix A is not of full rank in event ", iev
        continue

    # eta0: abs(eta) in [0.0,1.0]
    # eta1: abs(eta) in [1.0,1.5]
    # eta2: abs(eta) in [1.5,2.0]
    # eta3: abs(eta) in [2.0,2.5]

    #If any of the flavours of the jets are not 0 or 5, go to the next event. This is because sigmas were only fitted for Flavour = 0 or Flavour = 5
    if any(x not in [0,5] for x in ev.Jet_hadronFlavour):
        if print_discriminating_reasons:
            print "Jet w Flavour not in [0,5] in event " , str(iev)
        continue
        
    regions = np.zeros(ev.nJet)
    for jets in xrange(ev.nJet):
        if (np.absolute(ev.Jet_eta[jets]) > 0.0 and np.absolute(ev.Jet_eta[jets]) < 1.0):
            regions[jets] = 0
        elif (np.absolute(ev.Jet_eta[jets]) > 1.0 and np.absolute(ev.Jet_eta[jets]) < 1.5):
            regions[jets] = 1
        elif (np.absolute(ev.Jet_eta[jets]) > 1.5 and np.absolute(ev.Jet_eta[jets]) < 2.0):
            regions[jets] = 2
        elif (np.absolute(ev.Jet_eta[jets]) > 2.0 and np.absolute(ev.Jet_eta[jets]) < 2.5):
            regions[jets] = 3
        else: 
            regions[jets] = 99
    if any(x == 99 for x in regions):
        if print_discriminating_reasons:
            print "Jet w Eta not in any of the 4 regions in event ", str(iev), regions, jet_etas
        continue

    V = V_matrix(regions, ev.nJet, ev.Jet_pt, ev.Jet_eta, ev.Jet_hadronFlavour, RootFile)
    if np.linalg.matrix_rank(V) !=V.shape[0]:
        print "Matrix V is singular in event ", iev
        continue

    # R = [(- cos (V_phi) * V_pt),(-sin(V_phi)* V_pt)]
    # L is a 2 x nJet matrix, w. the first row containing the sines of the Jet_phis, the second row containing the cosines of the Jet_phis

    L, R = L_matrix_and_R_vector(ev.nJet, ev.Jet_pt, ev.V_pt, ev.V_eta, ev.V_phi, ev.V_mass, ev.Jet_phi, ev.Jet_mass, ev.Jet_eta)
        
    #Theta are the estimated values for the jet pts
    Theta = LagrangianSolver(A, L, V, R, ev.Jet_pt)

    Theta_before = Theta

    Lorentzvectors_before = []
    for i in xrange(Theta.shape[1]):
        v = ROOT.TLorentzVector()
        v.SetPtEtaPhiM(Theta[0,i],ev.Jet_eta[i], ev.Jet_phi[i], ev.Jet_mass[i])
        Lorentzvectors_before.append(v)

    lepton_vector = ROOT.TLorentzVector()
    lepton_vector.SetPtEtaPhiM(ev.V_pt, ev.V_eta, ev.V_phi, ev.V_mass)
    Lorentzvectors_before.append(lepton_vector)

    higgs_vector_before = []
    for i in ev.hJidx:
        v = ROOT.TLorentzVector()
        v.SetPtEtaPhiM(Theta[0,i], ev.Jet_eta[i], ev.Jet_phi[i], ev.Jet_mass[i])
        higgs_vector_before.append(v)
    
    addup = ROOT.TLorentzVector()
    for vector in higgs_vector_before:
        addup += vector
     
    contains_negativevalue = False

    if any(x < 0 for x in Theta[0,:].tolist()[0]):
        contains_negativevalue = True

#    if ev.nJet > 2:
#        if contains_negativevalue == True:
#            test = find_minimal_array(Theta[0,:].tolist()[0], ev.Jet_pt, V_matrix(regions, ev.nJet, ev.Jet_pt, ev.Jet_eta, ev.Jet_hadronFlavour, RootFile))
#            print test, "This is the Theta with minimal ChiSquare"
#            hue = raw_input()

    iteration_counter = 0.0
    positive_indices = []

    for i in xrange(Theta[0,:].shape[1]):
        if Theta[0,i] > 0:
            positive_indices.append(i)
    
    still_negative_value_left = []

    #Note that if contains_negativevalue isnt True on the first try, then the loop won't start at all
    while contains_negativevalue == True:
        
        contains_negativevalue = False

        for i in xrange(Theta[0,:].shape[1]):
            if Theta[0,i] < 0:
                contains_negativevalue = True

        positive_indices_scope = []
        for i in xrange(Theta[0,:].shape[1]):
            if Theta[0,i] > 0:
                positive_indices_scope.append(i)
                
        positive_indices = positive_indices_scope

        cur_V = V_matrix(regions, ev.nJet, ev.Jet_pt, ev.Jet_eta, ev.Jet_hadronFlavour, RootFile)
        if np.linalg.matrix_rank(cur_V) != cur_V.shape[0]:
            break_or_not = 1.0
            break

        positive_indices = find_minimal_array(Theta[0,:].tolist()[0], ev.Jet_pt, cur_V) 

#        new_theta_scope = [ev.Jet_pt[idx] for idx in positive_indices]
        new_theta_scope = [Theta[0,idx] for idx in positive_indices]
        new_regions_scope = [regions[idx] for idx in positive_indices]
        new_etas_scope = [ev.Jet_eta[idx] for idx in positive_indices]
        new_flavours_scope = [ev.Jet_hadronFlavour[idx] for idx in positive_indices]
        new_masses_scope = [ev.Jet_mass[idx] for idx in positive_indices]
        new_phis_scope = [ev.Jet_phi[idx] for idx in positive_indices] 
        


        break_or_not = 0.0
        if len(positive_indices) < 2:
#            print "--------------------"
#            print "Iterative process deleted too many jets, only ", len(positive_indices), " left in event ", iev, " after ", iteration_counter, " iterations"
            break_or_not = 1.0
        
        if break_or_not == 1.0:
            break

        A = A_matrix(a,b, len(positive_indices))
        V = V_matrix(new_regions_scope, len(positive_indices), new_theta_scope, new_etas_scope, new_flavours_scope, RootFile)
        if np.linalg.matrix_rank(V) != V.shape[0]:
            break_or_not == 1.0
            break
        L,R = L_matrix_and_R_vector(len(positive_indices), new_theta_scope, ev.V_pt, ev.V_eta, ev.V_phi, ev.V_mass, new_phis_scope, new_masses_scope, new_etas_scope)
        Theta = LagrangianSolver(A,L,V,R, new_theta_scope)
        
        iteration_counter += 1
    
#    print "There were ", iteration_counter, " iterations."

    for pt in Theta[0,:].tolist()[0]:
        if pt < 0:
            still_negative_value_left = True
    if still_negative_value_left:
        continue
        
    Lorentzvectors = []
    for i in xrange(len(positive_indices)):
        v = ROOT.TLorentzVector()
        v.SetPtEtaPhiM(Theta[0,i], [ev.Jet_eta[idx] for idx in positive_indices][i], [ev.Jet_phi[idx] for idx in positive_indices][i], [ev.Jet_mass[idx] for idx in positive_indices][i])
        Lorentzvectors.append(v)


    lepton_vector = ROOT.TLorentzVector()
    lepton_vector.SetPtEtaPhiM(ev.V_pt, ev.V_eta, ev.V_phi, ev.V_mass)
    Lorentzvectors.append(lepton_vector)

    px = 0.0
    py = 0.0

    for vector in Lorentzvectors_before:
        px += vector.Px()
        py += vector.Py()
#    print " "
#    print "Is pT conservation satisfied before iteration?"
#    print " "
#    print "px_before = ", px
#    print "py_before = ", py
#    print " "
    print Theta_before, "Theta before iteration"
    print Theta, "Theta after iteration"
    
    px = 0.0
    py = 0.0

    for i in xrange(len(Lorentzvectors)):
        px += Lorentzvectors[i].Px()
        py += Lorentzvectors[i].Py()
#    print " "
#    print "Fitted jets in event ", iev
#    print "Is pT conservation satisfied after iteration process?"
#    print " "
#    print "px_afterwards = ", px
#    print "py_afterwards = ", py
#    print " "

    higgs_vector_m = ROOT.TLorentzVector()
    higgs_vector = ROOT.TLorentzVector()
    hJidx_still_there = []
    for i in xrange(len(ev.hJidx)):
        v = ROOT.TLorentzVector()

        for j in xrange(len(ev.hJidx)):
            if ev.hJidx[i] in positive_indices:
                hJidx_still_there.append(True)
            else: hJidx_still_there.append(False)
        if all(hJidx_still_there):
            v.SetPtEtaPhiM(Theta[0,i], ev.Jet_eta[ev.hJidx[i]], ev.Jet_phi[ev.hJidx[i]], ev.Jet_mass[ev.hJidx[i]])
            higgs_vector += v
    
#        else: 
#            print "Deleted a hJidx jet"
        v.SetPtEtaPhiM(ev.Jet_pt[ev.hJidx[i]], ev.Jet_eta[ev.hJidx[i]], ev.Jet_phi[ev.hJidx[i]], ev.Jet_mass[ev.hJidx[i]])
        higgs_vector_m += v

    if any(x==False for x in hJidx_still_there):
        continue
            
    print addup.M(), "Estimated Higgs mass after first estimation"
    print "Estimated Higgs mass: ", higgs_vector.M()
    print "Measured Higgs mass: ", higgs_vector_m.M()
    print "---"

    no_fitted_events += 1

    phi_display = []
    pt_display = []
    eta_display = []
    for j in xrange(ev.nJet):
        phi_display.append(ev.Jet_phi[j])
        pt_display.append(ev.Jet_pt[j])
        eta_display.append(ev.Jet_eta[j])

    nJet[0] = ev.nJet
    est_mass[0] = higgs_vector.M()
    mea_mass[0] = higgs_vector_m.M()

    for jet in xrange(len(positive_indices)):
        estimates[jet] = Theta[0,jet]
        measurements[jet] = ev.Jet_pt[jet]
        etas[jet] = ev.Jet_eta[jet]        
        if jet in ev.hJidx:
            higgs_tag[jet] = 1
        else:
            higgs_tag[jet] = 0
#    print Theta[0,:], "Theta after Iteration"
#    print pt_display[0:ev.nJet]

    tree.Fill()
#    print "----------------------------------------------------"
print "We found ",no_fitted_events, " fitting events"    

out.cd()
tree.Write("", ROOT.TObject.kOverwrite)
out.Close()
