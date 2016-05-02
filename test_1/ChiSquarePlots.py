import numpy as np
import ROOT

ROOT.gStyle.SetOptStat(0)
RootFile = ROOT.TFile("ChiSquareFits_NeutrinoAdded_Full_only_Flavour5_hJCidx.root")
Tree = ROOT.gDirectory.Get('ChiSquareFits')

for no_jets in [2,3,4]:

    histo_1 = Tree.Draw("est_mass>>h1(100, 0, 250)", "nJet ==" + str(no_jets))
    histo_2 = Tree.Draw("mea_mass>>h2(100, 0, 250)", "nJet ==" + str(no_jets))
    histo_3 = Tree.Draw("(est_mass*(est_mass > 0) + 2*(est_mass == 0)*mea_mass + (est_mass > 0)*mea_mass)/2>>h3(100, 0, 250)", "nJet ==" + str(no_jets))

    h1 = ROOT.gDirectory.Get("h1")
    h2 = ROOT.gDirectory.Get("h2")
    h3 = ROOT.gDirectory.Get("h3")

    h1.SetTitle("nJet == " + str(no_jets))
    h2.SetTitle("nJet == " + str(no_jets))
    h3.SetTitle("nJet == " + str(no_jets))

    #Make Overflow visible by adding its contents to the last bin of the histogram
    no_bins_1 = h1.GetNbinsX()
    no_bins_2 = h2.GetNbinsX()
    no_bins_3 = h3.GetNbinsX()

    h1.AddBinContent(no_bins_1, h1.GetBinContent(no_bins_1+1))
    h2.AddBinContent(no_bins_2, h2.GetBinContent(no_bins_2+1))
    h3.AddBinContent(no_bins_3, h3.GetBinContent(no_bins_3+1))    

    h1.SetLineColor(4)
    h2.SetLineColor(2)
    h3.SetLineColor(6)

    c = ROOT.TCanvas()
    c.SetTitle("nJet = " + str(no_jets))

    h3.Draw()
    h1.Draw("SAME")
    h2.Draw("SAME")

    leg = ROOT.TLegend(.55, .5, .9, .7)
    leg.AddEntry(h1, "Estimated mass")
    leg.AddEntry(h2, "Measured mass")
    leg.AddEntry(h3, "(est_mass + mea_mass)/2")
    
    leg.Draw()

    img = ROOT.TImage.Create()
    img.FromPad(c)
    img.WriteImage("no_jets_"+str(no_jets)+".png")

histo_4 = Tree.Draw("est_mass>>h1", "nJet > 4")
histo_5 = Tree.Draw("mea_mass>>h2", "nJet > 4", "SAME")
histo_6 = Tree.Draw("(est_mass*(est_mass > 0) + 2*(est_mass == 0)*mea_mass + (est_mass > 0)*mea_mass)/2>>h3", "est_mass != 0 && nJet > 4", "SAME")

h1 = ROOT.gDirectory.Get("h1")
h2 = ROOT.gDirectory.Get("h2")
h3 = ROOT.gDirectory.Get("h3")

h1.SetTitle("nJet > 4")
h2.SetTitle("nJet > 4")
h3.SetTitle("nJet > 4")

h1.SetLineColor(4)
h2.SetLineColor(2)
h3.SetLineColor(6)

c = ROOT.TCanvas()

h3.Draw()
h1.Draw("SAME")
h2.Draw("SAME")

leg = ROOT.TLegend(.55, .5, .9, .7)
leg.SetHeader("nJet > 4 jets")
leg.AddEntry(h1, "Estimated mass")
leg.AddEntry(h2, "Measured mass")
leg.AddEntry(h3, "(est_mass + mea_mass)/2")

leg.Draw()

img = ROOT.TImage.Create()
img.FromPad(c)
img.WriteImage("no_jets_greaterthan_4.png")

