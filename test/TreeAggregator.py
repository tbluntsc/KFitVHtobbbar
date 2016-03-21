import os
import ROOT
import numpy as n

#Note the _small in the name of the root file: to try out new features
#loading all 25 files takes too long, so a smaller version of the root
#file will be created, only reading one of the root file in "address"
out = ROOT.TFile("Egen_to_Ereco_small.root", "UPDATE")
address = "dcap://t3se01.psi.ch:22125////pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-higgs/run2/VHBBHeppyV20/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/VHBB_HEPPY_V20_ZH_HToBB_ZToLL_M125_13TeV_powheg_Py8__fall15MAv2-pu25ns15v1_76r2as_v12-v1/160209_172236/0000/"

#Create an empty Tree 
tree = ROOT.TTree("Egen_to_Ereco", "Tree to compare RECO to GEN energies")

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

tree.Branch('RecoPt', Jet_pt, 'GenPt/F')
tree.Branch('GenPt', Jet_mcPt, 'RecoPt/F')
tree.Branch('RecoEta', Jet_eta, 'GenEta/F')
tree.Branch('RecoPhi', Jet_phi, 'GenPhi/F')
tree.Branch('GenEta', Jet_mcEta, 'RecoEta/F')
tree.Branch('GenPhi', Jet_mcPhi, 'RecoPhi/F')
tree.Branch('Flavour',Jet_hadronFlavour, 'GenFlavour/I')
tree.Branch('Jet_pt_reg', Jet_pt_reg, 'Jet_pt_reg/F')
tree.Branch('Mass', Jet_mass, 'GenMass/F')
tree.Branch('nJet', nJet, 'nJet/I')

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

print "Total number of entries: ", chain.GetEntries()

#To speed up runtime 3e+3 instead of 1e+11
for iev in range(int( min(1e+11, chain.GetEntries()))):
    chain.GetEntry(iev)
    ev = chain
    if iev%1000 == 0:
        print "Processing event ", iev+1
    for jets in xrange(ev.nJet):

        Jet_pt[0] = ev.Jet_pt[jets]
        Jet_hadronFlavour[0] = ev.Jet_hadronFlavour[jets]
        Jet_mcPt[0] = ev.Jet_mcPt[jets]
        Jet_eta[0] = ev.Jet_eta[jets]
        Jet_phi[0] = ev.Jet_phi[jets]
        Jet_mass[0] = ev.Jet_mass[jets]
        Jet_mcEta[0] = ev.Jet_mcEta[jets]
        Jet_mcPhi[0] = ev.Jet_mcPhi[jets]
        Jet_pt_reg[0] = ev.Jet_pt_reg[jets]
        nJet[0] = ev.nJet
        tree.Fill()
#    Jet_pt = ev.Jet_pt
#   Jet_hadronFlavour[0] = ev.Jet_hadronFlavour
#    Jet_mcPt = ev.Jet_mcPt
#    Jet_eta = ev.Jet_eta
#    Jet_phi = ev.Jet_phi
#    Jet_mass = ev.Jet_mass
#    Jet_mcEta = ev.Jet_mcEta
#    Jet_mcPhi = ev.Jet_mcPhi
#    Jet_pt_reg = ev.Jet_pt_reg

#for i in xrange(100):
#    interval = [i*5.0, (i+1)*5.0]


out.cd()
tree.Write("", ROOT.TObject.kOverwrite)
out.Close()
