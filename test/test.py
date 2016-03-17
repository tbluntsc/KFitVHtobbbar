import os
import ROOT
import numpy as n

#t = f.Get("tree")

xsec = 6025.2/1.23
tree = ROOT.TTree("tree_Z_inclusive", "tree title")
weight = 1.0

#This automatically makes all histograms of the TH1-type save the sum of the square of errors as well
ROOT.TH1.SetDefaultSumw2(True)
#This creates a root file called BFilter.root. The option update 
out = ROOT.TFile("BFilter.root", "UPDATE")

#initialize a bunch of zeros as an array of length 1
weight_ = n.zeros(1, dtype=float)
ttCls_ = n.zeros(1, dtype=float)
nGenStatus2bHad_ = n.zeros(1, dtype=float)
lheNb_ = n.zeros(1, dtype=float)
lheVpt_ = n.zeros(1, dtype=float)

#syntax: name of Branch, address of first item of structure, leaflist: i.e. list of variable names concatenated w all the variable types (see documentation for nomenclature of variable types ( D stands for 64 bit floating point )
tree.Branch('weight', weight_, 'weight/D')
tree.Branch('ttCls', ttCls_, 'ttCls/D')
tree.Branch('nGenStatus2bHad', nGenStatus2bHad_, 'nGenStatus2bHad/D')
tree.Branch('lheNb', lheNb_, 'lheNb/D')
tree.Branch('lheV_pt', lheVpt_, 'lheV_pt/D')

chain = ROOT.TChain('tree')

f = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125///pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-higgs/run2/VHBBHeppyV14/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V14_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151025_092939/0000/tree_1.root")
CountWeighted = 0
f.cd()
count = f.Get("CountWeighted")
CountWeighted += count.GetBinContent(1)
f.Close()
chain.AddFile("dcap://t3se01.psi.ch:22125///pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-higgs/run2/VHBBHeppyV14/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V14_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151025_092939/0000/tree_1.root")

weight = 10./(CountWeighted/xsec)*1000
print "Total processed events: ", CountWeighted

chain.SetBranchStatus("*", False)
chain.SetBranchStatus("ttCls", True)
chain.SetBranchStatus("nGenStatus2bHad", True)
chain.SetBranchStatus("lheNb", True)
chain.SetBranchStatus("lheV_pt", True)

print chain.GetEntries()
for iev in range( min(1e+11, chain.GetEntries())):
    chain.GetEntry(iev)
    ev = chain
    if iev%100 == 0:
        print "Processing", iev
    weight_[0] = weight
    ttCls_[0] = ev.ttCls
    nGenStatus2bHad_[0] = ev.nGenStatus2bHad
    lheNb_[0] = ev.lheNb
    lheVpt_[0] = ev.lheV_pt
    tree.Fill()

out.cd()
tree.Write("", ROOT.TObject.kOverwrite)
out.Close()
