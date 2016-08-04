import os
import ROOT
import numpy as n

#Note the _small in the name of the root file: to try out new features
#loading all 25 files takes too long, so a smaller version of the root
#file will be created, only reading one of the root file in "address"
#out = ROOT.TFile("Egen_to_Ereco_NeutrinosAdded_full_onlyFlavour5.root", "UPDATE")
out = ROOT.TFile("Resolutions/FlattenedTree_NOnuPt.root", "UPDATE")
out1 = ROOT.TFile("Resolutions/FlattenedTree_vLeptons.root", "UPDATE")
#address = "dcap://t3se01.psi.ch:22125////pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-higgs/run2/VHBBHeppyV20/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/VHBB_HEPPY_V20_ZH_HToBB_ZToLL_M125_13TeV_powheg_Py8__fall15MAv2-pu25ns15v1_76r2as_v12-v1/160209_172236/0000/"
address = "root://188.184.38.46:1094//store/group/phys_higgs/hbb/ntuples/V21/user/arizzi/VHBBHeppyV21/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/VHBB_HEPPY_V21_ZH_HToBB_ZToLL_M125_13TeV_powheg_Py8__fall15MAv2-pu25ns15v1_76r2as_v12-v1/160316_150654/0000/"
#address = "root://188.184.38.46:1094//store/group/phys_higgs/hbb/ntuples/V21/user/arizzi/VHBBHeppyV21/WJetsToQQ_HT180_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V21_WJetsToQQ_HT180_13TeV-madgraphMLM-Py8__fall15MAv2-pu25ns15v1_76r2as_v12-v1/160321_091519/0000"

#Create an empty Tree 
tree = ROOT.TTree("Egen_to_Ereco", "Tree to compare RECO to GEN energies")
tree1 = ROOT.TTree("tree", "vLeptons_pt and vLeptons_mcPt")

#Initialize desired Branches
Jet_pt_reg = n.zeros(1, dtype = n.float32())
Jet_pt = n.zeros(1, dtype = n.float32())
Jet_eta = n.zeros(1, dtype = n.float32())
Jet_phi = n.zeros(1, dtype = n.float32())
Jet_mass = n.zeros(1, dtype = n.float32())
Jet_mcPt = n.zeros(1, dtype = n.float32())
Jet_mcEta = n.zeros(1, dtype = n.float32())
Jet_mcPhi = n.zeros(1, dtype = n.float32())
Jet_hadronFlavour = n.zeros(1, dtype = int)
nJet = n.zeros(1, dtype = int)
vLeptons_pt = n.zeros(1, dtype = n.float32())
vLeptons_mcPt = n.zeros(1, dtype = n.float32())
V_pt = n.zeros(1, dtype = n.float32())
V_eta = n.zeros(1, dtype = n.float32())
V_phi = n.zeros(1, dtype = n.float32())
Jet_pt1 = n.zeros(1, dtype = n.float32())

tree.Branch('Jet_pt_reg', Jet_pt_reg, 'Jet_pt_reg/F')
tree.Branch('Jet_pt', Jet_pt, 'Jet_pt/F')
tree.Branch('Jet_eta', Jet_eta, 'Jet_eta/F')
tree.Branch('Jet_phi', Jet_phi, 'Jet_phi/F')
tree.Branch('Jet_mass', Jet_mass, 'Jet_mass/F')
tree.Branch('Jet_mcPt', Jet_mcPt, 'Jet_mcPt/F')
tree.Branch('Jet_mcEta', Jet_mcEta, 'Jet_mcEta/F')
tree.Branch('Jet_mcPhi', Jet_mcPhi, 'Jet_mcPhi/F')
tree.Branch('Jet_hadronFlavour',Jet_hadronFlavour, 'Jet_hadronFlavour/I')
tree.Branch('nJet', nJet, 'nJet/I')

tree1.Branch('vLeptons_pt', vLeptons_pt, 'vLeptons_pt/F')
tree1.Branch('vLeptons_mcPt', vLeptons_mcPt, 'vLeptons_mcPt/F')
tree1.Branch('Jet_pt', Jet_pt1, 'Jet_pt/F')

#Create array w. strings corresponding to tree names (i.e. tree_1.root, tree_2.root etc.)
no_of_files = 30
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
chain.SetBranchStatus("Jet_pt", True)
chain.SetBranchStatus("Jet_eta", True)
chain.SetBranchStatus("Jet_phi", True)
chain.SetBranchStatus("Jet_mass", True)
chain.SetBranchStatus("Jet_mcPt", True)
chain.SetBranchStatus("Jet_mcEta", True)
chain.SetBranchStatus("Jet_mcPhi", True)
chain.SetBranchStatus("Jet_pt_reg", True)
chain.SetBranchStatus("nJet", True)
chain.SetBranchStatus("nGenJet", True)
chain.SetBranchStatus("Jet_mcIdx", True)
chain.SetBranchStatus("GenJet_wNuPt", True)
chain.SetBranchStatus("vLeptons_pt", True)
chain.SetBranchStatus("vLeptons_mcPt", True)
chain.SetBranchStatus("V_eta", True)
chain.SetBranchStatus("V_phi", True)
chain.SetBranchStatus("V_pt", True)

print "Total number of entries: ", chain.GetEntries()

#To speed up runtime 3e+3 instead of 1e+11
for iev in range(int( min(1e+11, chain.GetEntries()))):
    chain.GetEntry(iev)
    ev = chain
    if iev%1000 == 0:
        print "Processing event ", iev+1

    # Jet_mcPt does not contain any information about possible neutrinos that may have been emitted
    # GenJet_wNuPt however is just that. This loop fils up new_Jet_mcPt w. GenJet_wNuPt if additional information is available
    # If the index is out of bounds or if it is -1 Jet_mcPt will be assigned 0 (and thus ignored in GaussianFit.py etc)

    new_Jet_mcPt = [ev.Jet_mcPt[i] for i in xrange(ev.nJet)]

#    for i in xrange(ev.nJet):
#        if ev.Jet_hadronFlavour[i] != 0:
#            if ev.Jet_mcIdx[i] == -1 or ev.Jet_mcIdx[i] > (len(ev.GenJet_wNuPt)-1):
#                new_Jet_mcPt[i] = 0
#            else:
#                new_Jet_mcPt[i] = ev.GenJet_wNuPt[ev.Jet_mcIdx[i]]

    for jets in xrange(ev.nJet):

        Jet_pt[0] = ev.Jet_pt[jets]
        Jet_hadronFlavour[0] = ev.Jet_hadronFlavour[jets]
        Jet_mcPt[0] = new_Jet_mcPt[jets]
        Jet_eta[0] = ev.Jet_eta[jets]
        Jet_phi[0] = ev.Jet_phi[jets]
        Jet_mass[0] = ev.Jet_mass[jets]
        Jet_mcEta[0] = ev.Jet_mcEta[jets]
        Jet_mcPhi[0] = ev.Jet_mcPhi[jets]
        Jet_pt_reg[0] = ev.Jet_pt_reg[jets]
        nJet[0] = ev.nJet

        tree.Fill()

    for lepton in xrange(len(ev.vLeptons_pt)):
        vLeptons_pt[0] = ev.vLeptons_pt[lepton]
        vLeptons_mcPt[0] = ev.vLeptons_mcPt[lepton]

        tree1.Fill()

out.cd()
tree.Write("", ROOT.TObject.kOverwrite)
out.Close()

out1.cd()
tree1.Write("", ROOT.TObject.kOverwrite)
out1.Close()
