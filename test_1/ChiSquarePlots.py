import numpy as np
import ROOT

RootFile = ROOT.TFile("ChiSquareFits.root")
Tree = ROOT.gDirectory.Get('ChiSquareFits')

for no_jets in [2,3,4]:

    histo_1 = Tree.Draw("est_mass>>h1", "est_mass < 300 && nJet ==" + str(no_jets))
    histo_2 = Tree.Draw("mea_mass>>h2", "est_mass < 300 && nJet ==" + str(no_jets))
    histo_3 = Tree.Draw("(est_mass + mea_mass)/2>>h3", "est_mass < 300 && nJet ==" + str(no_jets))

    h1 = ROOT.gDirectory.Get("h1")
    h2 = ROOT.gDirectory.Get("h2")
    h3 = ROOT.gDirectory.Get("h3")

    h1.SetLineColor(4)
    h2.SetLineColor(2)
    h3.SetLineColor(6)
    
    c = ROOT.TCanvas()

    h1.Draw()
    h2.Draw("SAME")
    h3.Draw("SAME")
    img = ROOT.TImage.Create()
    img.FromPad(c)
    img.WriteImage("no_jets_"+str(no_jets)+".png")
    
histo_4 = Tree.Draw("est_mass>>h1", "est_mass < 300 && nJet > 4")
histo_5 = Tree.Draw("mea_mass>>h2", "est_mass < 300 && nJet > 4", "SAME")
histo_6 = Tree.Draw("(est_mass+mea_mass)/2>>h3", "est_mass < 300 && nJet > 4", "SAME")

h1 = ROOT.gDirectory.Get("h1")
h2 = ROOT.gDirectory.Get("h2")
h3 = ROOT.gDirectory.Get("h3")

h1.SetLineColor(4)
h2.SetLineColor(2)
h3.SetLineColor(6)

c = ROOT.TCanvas()
    
h1.Draw()
h2.Draw("SAME")
h3.Draw("SAME")
img = ROOT.TImage.Create()
img.FromPad(c)
img.WriteImage("no_jets_greaterthan_4.png")
