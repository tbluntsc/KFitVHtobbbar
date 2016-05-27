import ROOT
import numpy as np
from numpy.linalg import inv

print_discriminating_reasons = False

#Opening RootFile containing the Sigmas Fits
RootFile = ROOT.TFile("V21_SigmasFits_test.root")

#Later in the code we will need to create three matrices to solve the Lagrangian system, A (model matrix), V (covariance matrix) & L (constraints matrix)
#Since these matrices will have to be built in each iteration of the loop over events a function is defined for each matrix which can be called inside the loop

def cut_off_one_element(array, higgsind):
    result = [[],[]]
    no_higgs_jets = sum(higgsind)
    for i in xrange(len(array)):
        cur_array = []
        cur_higgs = []
        for j in xrange(len(array)):
            if i != j:
                cur_array.append(array[j])
                cur_higgs.append(higgsind[j])
        if no_higgs_jets == sum(cur_higgs):    
            result[0].append(cur_array)
            result[1].append(cur_higgs)
    return result

def all_subarrays(nJet, higgs_indices):
    result = [[],[]]
    if nJet > 9:
        safe_result = []
        for i in xrange(nJet):
            safe_result.append(i)

        result[0] = safe_result
        result[1] = higgs_indices
        return [result]

    iter_counter = 0
    pts_indices = np.zeros(nJet, dtype = int)
    for i in xrange(nJet):
        pts_indices[i] = i

    initial_arrays = cut_off_one_element(pts_indices, higgs_indices)
    result[0].append(initial_arrays[0])
    result[1].append(initial_arrays[1])

    while(nJet - iter_counter > (sum(higgs_indices)+1)):
        cur_higgs = []
        cur_array = []

        for i in xrange(len(result[0][iter_counter][:])):

            intermediate_result = cut_off_one_element(result[0][iter_counter][i], result[1][iter_counter][i])

            for j in xrange(len(intermediate_result[0])):
                cur_array.append(intermediate_result[0][j])
                cur_higgs.append(intermediate_result[1][j])

        result[0].append(cur_array)
        result[1].append(cur_higgs)
        iter_counter += 1

    end_result_jets = []
    end_result_higgs = []
    for i in xrange(len(result[0])):
        cur_length = len(result[0][i][:])
        for j in xrange(cur_length):
            end_result_jets.append(result[0][i][j])
            end_result_higgs.append(result[1][i][j])

    combined_end_result = []
    for i in xrange(len(end_result_jets)):
        a = [end_result_jets[i], end_result_higgs[i]]
        combined_end_result.append(a)

    unique_combined_result = []
    for i in xrange(len(combined_end_result)):    
        if combined_end_result[i] not in unique_combined_result:
            unique_combined_result.append(combined_end_result[i])
            
    return unique_combined_result

def create_cut_permutations(Theta, higgs_indices):
    no_higgs_jets = sum(higgs_indices)
    result = []
    for i in xrange(len(Theta)):
        idx = []
        cur_higgs = []
        for j in xrange(len(Theta)):
            if j != i:
                cur_higgs.append(higgs_indices[j])
                idx.append(j)
        if sum(cur_higgs) == no_higgs_jets:
            result.append(idx)
    if len(result) < 1:
        return [0]
    return result

def ChiSquare(V, estimated_pt, original_pt):
    V_diag = np.diag(V).tolist()
    res = 0.0
    for idx in xrange(len(estimated_pt)):
        res += ((estimated_pt[idx] - original_pt[idx])**2)/(V_diag[idx])
    return res

def find_minimal_array(Theta, original_pt, V, higgs_indices):
    
    idx_list = create_cut_permutations(np.arange(len(Theta)), higgs_indices)
    if idx_list == [0]:
        return [0]
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
        
        jet_sigmas[jet] = myfunc.Eval(jet_pts[jet])*jet_pts[jet]

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
out = ROOT.TFile("V21_ChiSquareFitsUpgrade.root", "UPDATE")
#address = "dcap://t3se01.psi.ch:22125////pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-higgs/run2/VHBBHeppyV20/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/VHBB_HEPPY_V20_ZH_HToBB_ZToLL_M125_13TeV_powheg_Py8__fall15MAv2-pu25ns15v1_76r2as_v12-v1/160209_172236/0000/"
address = "root://188.184.38.46:1094//store/group/phys_higgs/hbb/ntuples/V21/user/arizzi/VHBBHeppyV21/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/VHBB_HEPPY_V21_ZH_HToBB_ZToLL_M125_13TeV_powheg_Py8__fall15MAv2-pu25ns15v1_76r2as_v12-v1/160316_150654/0000/"

#Create an empty Tree 
tree = ROOT.TTree("ChiSquareFits", "Tree containing results from kinematic fit")

#Initialize branches of tree
maxJet=100
nJet = np.zeros(1, dtype = int)
estimates = np.zeros(maxJet, dtype = np.float64())
measurements = np.zeros(maxJet, dtype = np.float64())
etas = np.zeros(maxJet, dtype = np.float64())
higgs_tag = np.zeros(maxJet, dtype = np.float64())
est_mass = np.zeros(1, dtype = np.float64())
mea_mass = np.zeros(1, dtype = np.float64())
Chisquare = np.zeros(1, dtype = np.float64())
nJet_after_it = np.zeros(1, dtype = int)

tree.Branch('nJet_after_it', nJet_after_it, 'nJet_after_it/I')
tree.Branch('nJet', nJet, 'nJet/I')
tree.Branch('estimates', estimates, "estimates[nJet_after_it]/D")
tree.Branch('measurements', measurements, "measurements[nJet]/D")
tree.Branch('etas', etas, "etas[nJet]/D")
tree.Branch('higgs_tag', higgs_tag, "higgs_tag[nJet_after_it]/D")
tree.Branch('est_mass', est_mass, "est_mass/D")
tree.Branch('mea_mass', mea_mass, "mea_mass/D")
tree.Branch('Chisquare', Chisquare, "Chisquare/D")

#Create array w. strings corresponding to tree names (i.e. tree_1.root, tree_2.root etc.)
no_of_files = 25
file_names = []
for i in range(1, no_of_files+1):
    file_names.append("tree_"+str(i)+".root")

print "Added ", len(file_names), "files"

#Loop over trees to chain them together
chain = ROOT.TChain("tree")

for file_name in file_names:
    f = ROOT.TFile.Open(address + "/" + file_name)
    if f==None or f.IsZombie():
        print "Address does not contain a tree"
        continue
    chain.AddFile(address + "/" + file_name)

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
chain.SetBranchStatus("hJCidx", True)

print "Total number of entries: ", chain.GetEntries()

no_fitted_events = 0.0

#To speed up runtime 3e+3 instead of 1e+11, i.e. only load the first 3000 events.
for iev in range(int( min(1e+11, chain.GetEntries()))):
    chain.GetEntry(iev)
    ev = chain
    
    #We wanna use Jet_pt_reg for Higgs jets (which we assume to be B-flavoured) and regular Jet_pt for the rest. 
    custom_pts = ev.Jet_pt
    for idx in ev.hJCidx:
        custom_pts[idx] = ev.Jet_pt_reg[idx]
    
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

    #Discard all entries w Higgs-tagged Jets that have PT < 20 or |eta| > 2.4
    higgs_bools = []
    for H_jets in xrange(len(ev.hJidx)):
        if (custom_pts[ev.hJidx[H_jets]] < 20) or (ev.Jet_eta[ev.hJidx[H_jets]] > 2.4) or (ev.Jet_eta[ev.hJidx[H_jets]] < -2.4) or ev.Jet_hadronFlavour[ev.hJidx[H_jets]] != 5:
            higgs_bools.append(True)
    if any(higgs_bools):
        if print_discriminating_reasons:
            print "Higgs jet w pt < 20 or |eta| > 2.4 in event ", str(iev)
        continue

    #Discard events where only one jet was produced (Is this unphysical? Or what's the problem? LC1LT seems to always be singular in that case)
    if ev.nJet == 1:
        if print_discriminating_reasons:
            print "Only one jet in event ", iev
        continue

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

    #Assign flavours to the jets. Higgs jets-> 5, else: 0
    Flavours = np.zeros(ev.nJet)
    
    for idx in ev.hJCidx:
        Flavours[idx] = 5
    
#    Flavours = ev.Jet_hadronFlavour

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

    for jets in xrange(ev.nJet):
        if regions[jets] != 99:
            histo_string = "Mean_" + str(int(regions[jets])) + "_" + str(int(jet_flavours[jet]))
            SigmasFile.cd()
            current_histo = ROOT.gDirectory.Get(histo_string)
            myfunc = current_histo.GetFunction("mean_func")
            custom_pts[jet] = custom_pts[jet] - (1.0 - myfunc.Eval(custom_pts[jet]))*custom_pts[jet]

    if any(x == 99 for x in regions):
        if print_discriminating_reasons:
            print "Jet w Eta not in any of the 4 regions in event ", str(iev), regions, jet_etas
        continue

    if len(ev.hJCidx) < 2:
        continue

    higgs_indices = []
    for i in xrange(ev.nJet):
        if i in ev.hJCidx:
            higgs_indices.append(1)
        else:
            higgs_indices.append(0)

    V = V_matrix(regions, ev.nJet, custom_pts, ev.Jet_eta, Flavours, RootFile)
    if np.linalg.matrix_rank(V) !=V.shape[0]:
        print "Matrix V is singular in event ", iev
        continue

    # R = [(- cos (V_phi) * V_pt),(-sin(V_phi)* V_pt)]
    # L is a 2 x nJet matrix, w. the first row containing the sines of the Jet_phis, the second row containing the cosines of the Jet_phis

    L, R = L_matrix_and_R_vector(ev.nJet, custom_pts, ev.V_pt, ev.V_eta, ev.V_phi, ev.V_mass, ev.Jet_phi, ev.Jet_mass, ev.Jet_eta)
        
    #Theta are the estimated values for the jet pts
    Theta = LagrangianSolver(A, L, V, R, custom_pts)
    Theta_before = Theta
#    print "-----------"
#    print [ev.Jet_pt[i] for i in xrange(ev.nJet)], "Original PTs" 
#    print Theta, "First estimation and its ChiSquare: ",  ChiSquare(V, Theta[0,:].tolist()[0], custom_pts)

    all_arrays = all_subarrays(ev.nJet, higgs_indices)


    lowest_ChiSquare = ChiSquare(V, Theta[0,:].tolist()[0], custom_pts)
    max_prob = ROOT.TMath.Prob(lowest_ChiSquare, ev.nJet)
    corresponding_result = Theta
    corresponding_indices = np.zeros(ev.nJet, dtype = int)

    for i in xrange(ev.nJet):
        corresponding_indices[i] = i
    
    if all_arrays != []:
        for index_arrays in all_arrays:
        
            iteration_regions = [regions[i] for i in index_arrays[0]]
            iteration_nJet = len(iteration_regions)
            iteration_pts = [custom_pts[i] for i in index_arrays[0]]
            iteration_etas = [ev.Jet_eta[i] for i in index_arrays[0]]
            iteration_Flavours = [Flavours[i] for i in index_arrays[0]]
            iteration_phis = [ev.Jet_phi[i] for i in index_arrays[0]]
            iteration_masses = [ev.Jet_mass[i] for i in index_arrays[0]]
        
            V = V_matrix(iteration_regions, iteration_nJet, iteration_pts, iteration_etas, iteration_Flavours, RootFile)
            L, R = L_matrix_and_R_vector(iteration_nJet, iteration_pts, ev.V_pt, ev.V_eta, ev.V_phi, ev.V_mass, iteration_phis, iteration_masses, iteration_etas)
            A = A_matrix(a,b, len(index_arrays[0]))
        
            iteration_result = LagrangianSolver(A, L, V, R, iteration_pts)
#        print iteration_result," result ", i ," and its ChisQuare: ", ChiSquare(V, iteration_result[0,:].tolist()[0], iteration_pts)

            iteration_ChiSquare = ChiSquare(V, iteration_result[0,:].tolist()[0], iteration_pts)
            iteration_prob = ROOT.TMath.Prob(iteration_ChiSquare, iteration_nJet)

#        print iteration_prob, "Yo probability for arrays ", index_arrays
            neg_value_or_no = False

#        print "Iteration result: ", iteration_result[0,:].tolist()[0]

            for pt in iteration_result[0,:].tolist()[0]:
                if pt < 0:
                    neg_value_or_no = True

            if not neg_value_or_no:
                if iteration_prob > max_prob:
                    corresponding_result = iteration_result
                    corresponding_indices = index_arrays[0]
#    print ">>>>>>>>>>>>>>>>>>>>>>>>>>>Below, the result<<<<<<<<<<<<<<<<<<<<<<<<<"
#    print corresponding_result
#    print "and the original pts ", [custom_pts[i] for i in xrange(len(custom_pts))] 
 #   print corresponding_indices
 #   print "----------"
    
    new_etas = [ev.Jet_eta[i] for i in corresponding_indices]
    new_phis = [ev.Jet_phi[i] for i in corresponding_indices]

#    print corresponding_indices, "indices of result"

    scale_factors = np.zeros(len(corresponding_indices))
    for i in xrange(len(scale_factors)):
        scale_factors[i] = corresponding_result[0,i]/custom_pts[corresponding_indices[i]]

    new_masses = np.zeros(len(corresponding_indices))
    for i in xrange(len(new_masses)):
        new_masses[i] = ev.Jet_mass[corresponding_indices[i]]*scale_factors[i]

#    print [ev.Jet_mass[i] for i in corresponding_indices]
#    print new_masses
        

    Lorentzvectors_before = []
    for i in xrange(Theta_before.shape[1]):
        v = ROOT.TLorentzVector()
        v.SetPtEtaPhiM(Theta_before[0,i],ev.Jet_eta[i], ev.Jet_phi[i], ev.Jet_mass[i])
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

    still_negative_value_left = False
 
    for pt in corresponding_result[0,:].tolist()[0]:
        if pt < 0:
            still_negative_value_left = True
            
    corresponding_result = corresponding_result[0,:].tolist()[0]
    Lorentzvectors = []
    for i in xrange(len(corresponding_result)):
        v = ROOT.TLorentzVector()
        v.SetPtEtaPhiM(corresponding_result[i], new_etas[i],  new_phis[i], new_masses[i])
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
#    print [ev.Jet_pt[i] for i in xrange(ev.nJet)], "Original pts"
#    print Theta_before, "Theta before iteration"
#    print Theta, "Theta after iteration"
#    print " "
    px = 0.0
    py = 0.0

    for i in xrange(len(Lorentzvectors)):
        px += Lorentzvectors[i].Px()
        py += Lorentzvectors[i].Py()
#    print " "
#    print "Is pT conservation satisfied after iteration process?"
#    print " "
    if (px > 0.001 or py > 0.001) and still_negative_value_left == False :
        print "px_afterwards = ", px
        print "py_afterwards = ", py
        print still_negative_value_left
        print " "
        
    higgs_vector_m = ROOT.TLorentzVector()
    higgs_vector = ROOT.TLorentzVector()
    
    for idx in ev.hJidx:
        cur_v = ROOT.TLorentzVector()
        cur_v.SetPtEtaPhiM(custom_pts[idx], ev.Jet_eta[idx], ev.Jet_phi[idx], ev.Jet_mass[idx])
        higgs_vector_m += cur_v

    for i in xrange(len(corresponding_indices)):
        cur_v = ROOT.TLorentzVector()
        if higgs_indices[corresponding_indices[i]] == 1:
            cur_v.SetPtEtaPhiM(corresponding_result[i], new_etas[i], new_phis[i], new_masses[i])
            higgs_vector += cur_v

#    print addup.M(), "Estimated Higgs mass after first estimation"
#    print "Estimated Higgs mass: ", higgs_vector.M()
#    print "Measured Higgs mass: ", higgs_vector_m.M()
 #   print len(final_indices), "nJet after it"
 #   print final_indices, "final indices"
    no_fitted_events += 1

    #Fill up the output tree with the results and some metadata
    nJet[0] = ev.nJet
    nJet_after_it[0] = len(corresponding_indices)
    mea_mass[0] = higgs_vector_m.M()

    if still_negative_value_left:
        Chisquare[0] = 0.0
        est_mass[0] = 0.0
        for i in xrange(nJet_after_it):
            estimates[i] = 0.0
            higgs_tag[i] = higgs_indices[i]


    else:
        V = V_matrix([regions[i] for i in corresponding_indices], len(corresponding_result), corresponding_result, [ev.Jet_eta[i] for i in corresponding_indices], [Flavours[i] for i in corresponding_indices], RootFile)
        Chisquare[0] = ChiSquare(V, corresponding_result, [custom_pts[i] for i in corresponding_indices])
        est_mass[0] = higgs_vector.M()

        for i in xrange(nJet_after_it):
            estimates[i] = Theta[0,i]
            higgs_tag[i] = higgs_indices[i]
  
    for jet in xrange(ev.nJet):
        measurements[jet] = custom_pts[jet]
        etas[jet] = ev.Jet_eta[jet]  

#    print '--------------------------------------------------'
    tree.Fill()

print "We found ",no_fitted_events, " fitting events"    
out.cd()
tree.Write("", ROOT.TObject.kOverwrite)
out.Close()
