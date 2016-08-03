import numpy as np
import ROOT

RootFile = ROOT.TFile("Resolutions/FlattenedTree_WITHnuPt.root")
Tree = ROOT.gDirectory.Get('Egen_to_Ereco')
ROOT.gStyle.SetOptStat(0)
ROOT.TGaxis.SetMaxDigits(3)

for Jet_pt in [25,50,75,100,125,150,175]:


    histo_1 = Tree.Draw("Jet_mcEta-Jet_eta>>h1(100, -0.1, 0.1)", "Jet_pt > " + str(Jet_pt) + "&& Jet_pt < " + str(Jet_pt+ 25))
    h1 = ROOT.gDirectory.Get("h1")

    h1.SetTitle("")

    c = ROOT.TCanvas(str(Jet_pt), str(Jet_pt), 500,650)
    c.SetRightMargin(0.04)
    c.SetTopMargin(0.21)
    c.SetLeftMargin(0.13)
    c.SetBottomMargin(0.15)

    pt = ROOT.TPaveText(0.2, 0.40, 0.45, 0.77, "NDC")
    pt.SetTextSize(0.021)
    pt.AddText("Difference between")
    pt.AddText("Reco- & Gen-level #eta")
    pt.AddLine(.0,.6,1.,.6)
    pt.AddText("POWHEG, PYTHIA8, MC")
    pt.AddText("ZH,H->bb,Z->ll")
    pt.AddText("|#vec{p_{T}}| #in ["+str(Jet_pt)+","+str(Jet_pt+25)+"]")

    x_axis = h1.GetXaxis()
    x_axis.SetTitle("GenEta - RecoEta (Radians)")
    x_axis.SetTitleSize(0.078)
    x_axis.SetTitleOffset(0.6)
    x_axis.SetLabelSize(0.028)

    y_axis = h1.GetYaxis()
    y_axis.SetTitle("Events")
    y_axis.SetTitleSize(0.078)
    y_axis.SetTitleOffset(0.6)
    y_axis.SetLabelSize(0.028)


    first = h1.FindFirstBinAbove(h1.GetMaximum()/2)
    last = h1.FindLastBinAbove(h1.GetMaximum()/2)
    fwhm = (h1.GetBinCenter(last) - h1.GetBinCenter(first))

    leg = ROOT.TLegend(0.13, 0.82, 0.94, 0.98)
    leg.SetTextSize(0.08)
    leg.AddEntry(h1, "Mean = " + str(round(h1.GetMean(),6)))
    leg.AddEntry(h1, "FWHM = " + str(fwhm))

    h1.Draw()
    leg.Draw()
    pt.Draw()

    img = ROOT.TImage.Create()
    img.FromPad(c)
    img.WriteImage("Eta_resolutions" + str(Jet_pt) + ".png")

    c.Close()
