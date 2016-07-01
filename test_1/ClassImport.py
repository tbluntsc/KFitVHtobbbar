import EstimationClass as Hue
import ROOT
import numpy as np
from numpy.linalg import inv

est = Hue.Estimation()

address = "root://188.184.38.46:1094//store/group/phys_higgs/hbb/ntuples/V21/user/arizzi/VHBBHeppyV21/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/VHBB_HEPPY_V21_ZH_HToBB_ZToLL_M125_13TeV_powheg_Py8__fall15MAv2-pu25ns15v1_76r2as_v12-v1/160316_150654/0000/"

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
difference_counter = 0.0

out = ROOT.TFile("Current_Root_Files/ClassTest.root", "UPDATE")
tree = ROOT.TTree("ChiSquareFits", "Tree containing results from kinematic fit")
nJet = np.zeros(1, dtype = int)
est_mass = np.zeros(1, dtype = np.float64())
mea_mass = np.zeros(1, dtype = np.float64())
mea_mass_reg = np.zeros(1, dtype = np.float64())

tree.Branch('nJet', nJet, 'nJet/I')
tree.Branch('est_mass', est_mass, "est_mass/D")
tree.Branch('mea_mass', mea_mass, "mea_mass/D")
tree.Branch('mea_mass_reg', mea_mass_reg, "mea_mass_reg/D")

no_discriminated_events = 0
no_fits = 0
no_unsuccessful_estimations = 0
#To speed up runtime 3e+3 instead of 1e+11, i.e. only load the first 3000 events.
for iev in range(int( min(1e+11, chain.GetEntries()))):
    chain.GetEntry(iev)
    if iev%5000 == 0:
        print "Sup dawg, processing event ", iev
#        print "and we have found ", no_fits, " fitting events so far"
#        print "and have discriminated ", no_discriminated_events, " events"
#        print "and have had ", no_unsuccessful_estimations, " unsuccesful estimations
#        print est.discr_lessthan2Higgsjets, "< 2 higgs jets"
#        print est.discr_Vtype, "Vtype"
#        print est.discr_Vpt , " Vpt"
#        print est.discr_lessthan2jets, "< 2 jets"
#        print est.discr_noresolution, "no resolution"
#        print est.discr_HiggsPtEtaFlavour, "HiggsPtEtaFlavour" 
    ev = chain

    est.Discriminate_and_custompts(ev)

    if est.Discriminate == False:
        est.mass(ev)
#        print est.estimation_successful, "<- whether estimation was successful or nah"
#        print est.estimated_mass, "estmass"
#        print est.measured_mass, "meamass"
#        print est.measured_mass_reg, "mea_mass_reg"
#        print "----------------------------------"

        if est.estimation_successful:
            est_mass[0] = est.estimated_mass
            mea_mass[0] = est.measured_mass
            mea_mass_reg[0] = est.measured_mass_reg
            nJet[0] = ev.nJet
            tree.Fill()
            no_fits += 1
        else: 
            no_unsuccessful_estimations += 1
    else:
        no_discriminated_events += 1

        
print "and we have found ", no_fits, " fitting events so far"
print "and have discriminated ", no_discriminated_events, " events"
print "and have had ", no_unsuccessful_estimations, " unsuccesful estimations"
print est.discr_lessthan2Higgsjets, "< 2 higgs jets"
print est.discr_Vtype, "Vtype"
print est.discr_Vpt , " Vpt"
print est.discr_lessthan2jets, "< 2 jets"
print est.discr_noresolution, "no resolution"
print est.discr_HiggsPtEtaFlavour, "HiggsPtEtaFlavour" 



out.cd()
tree.Write("", ROOT.TObject.kOverwrite)
out.Close()
