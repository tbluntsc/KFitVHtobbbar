import ROOT
import numpy as np
from numpy.linalg import inv

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
    if nJet > 8:
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

def ChiSquare(V, estimated_pt, original_pt):
    V_diag = np.diag(V).tolist()
    res = 0.0
    for idx in xrange(len(estimated_pt)):
        res += ((estimated_pt[idx] - original_pt[idx])**2)/(V_diag[idx])
    return res

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

        cur_sigma = myfunc.Eval(jet_pts[jet])
        if (cur_sigma < 0.00000001) & (cur_sigma > -0.00000001):
            cur_sigma = 0.001

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

def Estimated_pts(event, Resolutions, print_discriminating_reasons = False):

    did_we_estimate = False

    no_estimate_indices = np.zeros(event.nJet, dtype = int)
    for i in xrange(event.nJet):
        no_estimate_indices[i] = i

    #There can't be less than 2 Higgs jets
    if len(event.hJCidx) < 2:
        return event.Jet_pt, no_estimate_indices, did_we_estimate

    #Discard all entries w Vtype not equal to either 0 or 1. (1: V -> e+e- / e-v_e, 0: V -> mumu / mu v_mu )
    if (event.Vtype != 0) & (event.Vtype != 1):
        if print_discriminating_reasons:
            print "VType neither 0 or 1, Vtype = ", str(event.Vtype) 
        return event.Jet_pt, no_estimate_indices, did_we_estimate

    #Discard all entries w Vectorboson Momentum smaller than 50 Gevent
    if event.V_pt < 50:
        if print_discriminating_reasons:
            print "V_pt < 50, V_pt = ", str(event.V_pt)
        return event.Jet_pt, no_estimate_indices, did_we_estimate

    #Discard events where only one jet was produced
    if event.nJet == 1:
        if print_discriminating_reasons:
            print "Only one jet"
        return event.Jet_pt, no_estimate_indices, did_we_estimate

    # eta0: abs(eta) in [0.0,1.0]
    # eta1: abs(eta) in [1.0,1.5]
    # eta2: abs(eta) in [1.5,2.0]
    # eta3: abs(eta) in [2.0,2.5]

    regions = np.zeros(event.nJet)
    for jets in xrange(event.nJet):
        if (np.absolute(event.Jet_eta[jets]) > 0.0 and np.absolute(event.Jet_eta[jets]) < 1.0):
            regions[jets] = 0
        elif (np.absolute(event.Jet_eta[jets]) > 1.0 and np.absolute(event.Jet_eta[jets]) < 1.5):
            regions[jets] = 1
        elif (np.absolute(event.Jet_eta[jets]) > 1.5 and np.absolute(event.Jet_eta[jets]) < 2.0):
            regions[jets] = 2
        elif (np.absolute(event.Jet_eta[jets]) > 2.0 and np.absolute(event.Jet_eta[jets]) < 2.5):
            regions[jets] = 3
        else: 
            regions[jets] = 99

    if any(x == 99 for x in regions):
        if print_discriminating_reasons:
            print "Jet w Eta not in any of the 4 regions"
        return event.Jet_pt, no_estimate_indices, did_we_estimate
    
    #We wanna use Jet_pt_reg for Higgs jets (which we assume to be B-flavoured) and regular Jet_pt for the rest. 
    custom_pts = np.zeros(ev.nJet)

    for i in xrange(ev.nJet):
        custom_pts[i] = ev.Jet_pt[i]
    for idx in ev.hJCidx:
        custom_pts[idx] = ev.Jet_pt_reg[idx]


    #Assign flavours to the jets. Higgs jets-> 5, else: 0
    Flavours = np.zeros(event.nJet)
    for idx in event.hJCidx:
        Flavours[idx] = 5

    #Normalize mean of resolutions to 1
    for jets in xrange(event.nJet):

        histo_string = "Mean_" + str(int(regions[jets])) + "_" + str(int(Flavours[jets]))
        Resolutions.cd()
        current_histo = ROOT.gDirectory.Get(histo_string)
        myfunc = current_histo.GetFunction("mean_func")
        custom_pts[jets] = custom_pts[jets]/(myfunc.Eval(custom_pts[jets]))

    #Discard all entries w Higgs-tagged Jets that have PT < 20 or |eta| > 2.4
    higgs_bools = []
    for H_jets in xrange(len(event.hJCidx)):
        if (custom_pts[event.hJCidx[H_jets]] < 20) or (event.Jet_eta[event.hJCidx[H_jets]] > 2.4) or (event.Jet_eta[event.hJCidx[H_jets]] < -2.4):
            higgs_bools.append(True)
    if any(higgs_bools):
        if print_discriminating_reasons:
            print "Higgs jet w pt < 20 or |eta| > 2.4 in event ", str(ievent)
        return event.Jet_pt, no_estimate_indices, did_we_estimate


    higgs_indices = []
    for i in xrange(event.nJet):
        if i in event.hJCidx:
            higgs_indices.append(1)
        else:
            higgs_indices.append(0)

    ##### Building matrices & solving Lagrangian system #####

    # R = [(- cos (V_phi) * V_pt),(-sin(V_phi)* V_pt)]
    # L is a 2 x nJet matrix, w. the first row containing the sines of the Jet_phis, the second row containing the cosines of the Jet_phis
    #Theta are the estimated values for the jet pts

    a = 1.0
    b = 0.0

    A = A_matrix(a,b, event.nJet)

    if np.linalg.matrix_rank(A) != A.shape[0]:
        print "Matrix A is singular"
        return event.Jet_pt, no_estimate_indices, did_we_estimate

    V = V_matrix(regions, event.nJet, custom_pts, event.Jet_eta, Flavours, Resolutions)

    if np.linalg.matrix_rank(V) !=V.shape[0]:
        print "Matrix V is singular"
        return event.Jet_pt, no_estimate_indices, did_we_estimate
    L, R = L_matrix_and_R_vector(event.nJet, custom_pts, event.V_pt, event.V_eta, event.V_phi, event.V_mass, event.Jet_phi, event.Jet_mass, event.Jet_eta)


    Theta = LagrangianSolver(A, L, V, R, custom_pts)

    #List containing all possible index combinations w nJet, nJet-1, ..., len(ev.hJCidx) jets without deleting any Higgs jets.
    all_arrays = all_subarrays(event.nJet, higgs_indices)

    #Inialize first guess for best estimation
    lowest_ChiSquare = ChiSquare(V, Theta[0,:].tolist()[0], custom_pts)
    max_prob = ROOT.TMath.Prob(lowest_ChiSquare, event.nJet)
    corresponding_result = Theta
    corresponding_indices = np.zeros(event.nJet, dtype = int)
    for i in xrange(event.nJet):
        corresponding_indices[i] = i
    

    if all_arrays != []:
        for index_arrays in all_arrays:
            iteration_regions = [regions[i] for i in index_arrays[0]]
            iteration_nJet = len(iteration_regions)
            iteration_pts = [custom_pts[i] for i in index_arrays[0]]
            iteration_etas = [event.Jet_eta[i] for i in index_arrays[0]]
            iteration_Flavours = [Flavours[i] for i in index_arrays[0]]
            iteration_phis = [event.Jet_phi[i] for i in index_arrays[0]]
            iteration_masses = [event.Jet_mass[i] for i in index_arrays[0]]
        
            V = V_matrix(iteration_regions, iteration_nJet, iteration_pts, iteration_etas, iteration_Flavours, Resolutions)
            if np.linalg.matrix_rank(V) !=V.shape[0]:
                continue
            L, R = L_matrix_and_R_vector(iteration_nJet, iteration_pts, event.V_pt, event.V_eta, event.V_phi, event.V_mass, iteration_phis, iteration_masses, iteration_etas)
            if np.linalg.matrix_rank(L) !=L.shape[0]:
                continue
            A = A_matrix(a,b, len(index_arrays[0]))
            if np.linalg.matrix_rank(A) !=A.shape[0]:
                continue
        
            iteration_result = LagrangianSolver(A, L, V, R, iteration_pts)
            iteration_ChiSquare = ChiSquare(V, iteration_result[0,:].tolist()[0], iteration_pts)
            iteration_prob = ROOT.TMath.Prob(iteration_ChiSquare, iteration_nJet)

            neg_value_or_no = False
            for pt in iteration_result[0,:].tolist()[0]:
                if pt < 0:
                    neg_value_or_no = True

            if not neg_value_or_no:
                if iteration_prob > max_prob:
                    corresponding_result = iteration_result
                    corresponding_indices = index_arrays[0]
 
    new_etas = [event.Jet_eta[i] for i in corresponding_indices]
    new_phis = [event.Jet_phi[i] for i in corresponding_indices]

    still_negative_value_left = False
 
    for pt in corresponding_result[0,:].tolist()[0]:
        if pt < 0:
            still_negative_value_left = True

    if still_negative_value_left:
        return event.Jet_pt, no_estimate_indices, did_we_estimate
            
    corresponding_result = corresponding_result[0,:].tolist()[0]

#    print "Estimation succesful"
    did_we_estimate = True
    return corresponding_result, corresponding_indices, did_we_estimate


RootFile = ROOT.TFile("V21_SigmasFits_fixedNames.root")
address = "root://188.184.38.46:1094//store/group/phys_higgs/hbb/ntuples/V21/user/arizzi/VHBBHeppyV21/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/VHBB_HEPPY_V21_ZH_HToBB_ZToLL_M125_13TeV_powheg_Py8__fall15MAv2-pu25ns15v1_76r2as_v12-v1/160316_150654/0000/"
#Create array w. strings corresponding to tree names (i.e. tree_1.root, tree_2.root etc.)
no_of_files = 30
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

out = ROOT.TFile("Current_Root_Files/FunctionTest.root", "UPDATE")
tree = ROOT.TTree("ChiSquareFits", "Tree containing results from kinematic fit")
nJet = np.zeros(1, dtype = int)
est_mass = np.zeros(1, dtype = np.float64())
mea_mass = np.zeros(1, dtype = np.float64())

tree.Branch('nJet', nJet, 'nJet/I')
tree.Branch('est_mass', est_mass, "est_mass/D")
tree.Branch('mea_mass', mea_mass, "mea_mass/D")

#To speed up runtime 3e+3 instead of 1e+11, i.e. only load the first 3000 events.
for iev in range(int( min(1e+11, chain.GetEntries()))):
    chain.GetEntry(iev)
    ev = chain
    if iev%1000 == 0:
        print "Processing event ", iev

    pts, indices, save = Estimated_pts(ev, RootFile, print_discriminating_reasons = False)

    if not save:
        continue

    higgs_indices = np.zeros(ev.nJet, dtype = int)
    for idx in ev.hJCidx:
        higgs_indices[idx] = 1

    higgs_post_est = [higgs_indices[i] for i in indices]

    custom_pts = np.zeros(ev.nJet)
    for i in xrange(ev.nJet):
        if higgs_indices[i] == 1:
            custom_pts[i] = ev.Jet_pt_reg[i]
        else:
            custom_pts[i] = ev.Jet_pt[i]

    Flavours = np.zeros(ev.nJet, dtype= int)
    for idx in ev.hJCidx:
        Flavours[idx] = 5

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

        histo_string = "Mean_" + str(int(regions[jets])) + "_" + str(int(Flavours[jets])\
)
        RootFile.cd()
        current_histo = ROOT.gDirectory.Get(histo_string)
        myfunc = current_histo.GetFunction("mean_func")
        custom_pts[jets] = custom_pts[jets]/(myfunc.Eval(custom_pts[jets]))

    est_vector = ROOT.TLorentzVector()
    mea_vector = ROOT.TLorentzVector()

    scale_factors = np.zeros(len(indices))

    for i in xrange(len(scale_factors)):
        scale_factors[i] = pts[i]/custom_pts[indices[i]]

    new_masses = np.zeros(len(indices))
    for i in xrange(len(new_masses)):
        new_masses[i] = ev.Jet_mass[indices[i]]*scale_factors[i]

    for i in xrange(len(indices)):
        if higgs_post_est[i] == 1:
            v = ROOT.TLorentzVector()
            v.SetPtEtaPhiM(pts[i], ev.Jet_eta[indices[i]], ev.Jet_phi[indices[i]], new_masses[i])
            est_vector += v

    higgs_vector_m_reg = ROOT.TLorentzVector()
    higgs_vector_m = ROOT.TLorentzVector()

    for idx in ev.hJCidx:
        cur_v = ROOT.TLorentzVector()
        cur_v.SetPtEtaPhiM(ev.Jet_pt_reg[idx], ev.Jet_eta[idx], ev.Jet_phi[idx], ev.Jet_mass[idx])
        higgs_vector_m_reg += cur_v

        cur_v = ROOT.TLorentzVector()
        cur_v.SetPtEtaPhiM(ev.Jet_pt[idx], ev.Jet_eta[idx], ev.Jet_phi[idx], ev.Jet_mass[idx])
        higgs_vector_m += cur_v

    mea_mass[0] = higgs_vector_m.M()
    mea_mass_reg = higgs_vector_m_reg.M()
    nJet[0] = ev.nJet
    est_mass[0] = est_vector.M()
    tree.Fill()

out.cd()
tree.Write("", ROOT.TObject.kOverwrite)
out.Close()
