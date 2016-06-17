import ROOT
import numpy as np
from matplotlib import pyplot as plt
ROOT.gSystem.Load("libRooFit")

#Chisquare_values = ["0.1", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "4.0", "4.5", "5.0"]
Chisquare_values = ["0.01","0.05","0.1"]

differences = []
rel_uncertainties = []
errors = []
labeler = 0

for Chisquare in Chisquare_values:

    RootHistogram = ROOT.RooDataHist()

    c1 = ROOT.TCanvas("c1", "c1", 800,800)
    leg = ROOT.TLegend(0.64, 0.65, 0.9, 0.9)

    x = ROOT.RooRealVar("x-axis", "x", 20,180)
    Xp = ROOT.RooRealVar("Xp_sgn", "Xp",100.0 , 200.0)
    sP = ROOT.RooRealVar("sP_sgn", "sP", 0.0 , 50.0)
    xi = ROOT.RooRealVar("xi_sgn", "xi", -5.0, 5.0)
    rho1 = ROOT.RooRealVar("rho1_sgn", "rho1", -2.0, 2.0)
    rho2 = ROOT.RooRealVar("rho2_sgn", "rho2", -2.0, 2.0)
    RooBukinPdf = ROOT.RooBukinPdf("buk_pdf_sgn","RooBukinPdf for signal", x, Xp, sP, xi, rho1, rho2)
        
    frame = x.frame()
    frame.SetName("frame")
    frame.SetTitle("Higgs mass estimate w. a Bukin function fitted to it for Chisquare < " +Chisquare)

    for variable in ["est_mass", "mea_mass_reg"]:

        data = ROOT.TFile('Current_Root_Files/V21_ChiSquareFitsUpgrade_hJCidx.root')
        Tree = ROOT.gDirectory.Get('ChiSquareFits')

        #Get number of events before restricting Chisquare
        Tree.Draw(variable+">>h1(100,20,180)", "est_mass < 250 && est_mass != 0")
        h1 = ROOT.gDirectory.Get("h1")
 
        Roo_histo_before = ROOT.RooDataHist("higgs mass w/o Chisquare restriction", "higgs mass w/o Chisquare restriction", ROOT.RooArgList(x), h1, 1.0)
        counts_before = h1.Integral()

#        RooBukinPdf.fitTo(Roo_histo_before)        

        name = ROOT.RooFit.Name(Chisquare)
        
        if variable == "est_mass":
            color = ROOT.RooFit.LineColor(4)
            style_histo = ROOT.RooFit.MarkerStyle(1)
        if variable == "mea_mass_reg":
            color = ROOT.RooFit.LineColor(5)
            style_histo = ROOT.RooFit.MarkerStyle(2)

#        Roo_histo_before.plotOn(frame, color, name)
#        RooBukinPdf.plotOn(frame, color, name)

        cur_mean = str(Xp.getVal())
        cur_sig = str(sP.getVal())
        cur_unc = str(sP.getVal() / Xp.getVal())

#        leg.AddEntry(frame.getCurve(variable), variable)
#        leg.AddEntry(frame.getCurve(variable), "mean = " + cur_mean)
#        leg.AddEntry(frame.getCurve(variable), "sigma = " + cur_sig)
#        leg.AddEntry(frame.getCurve(variable), "rel-uncert. = " + cur_unc)

        #Create histogram w. restricted Chisquare and get the difference of events 
        Tree.Draw(variable+">>h2(100,20,180)", "est_mass < 250 && est_mass != 0  && Chisquare <" +Chisquare)
        h2 = ROOT.gDirectory.Get("h2")

        Roo_histo_1 = ROOT.RooDataHist("higgs mass w. Chisquare restriction ", "estimated higgs mass",  ROOT.RooArgList(x), h2, 1.0)
        counts_after = h2.Integral()
        differences.append(1 - counts_after/counts_before)

        RooBukinPdf.fitTo(Roo_histo_1)        
        if variable == "est_mass":
            rel_uncertainties.append(sP.getVal()/Xp.getVal())

        name = ROOT.RooFit.Name(Chisquare)
        
        if variable == "est_mass":
            color = ROOT.RooFit.LineColor(2)
            style_histo = ROOT.RooFit.MarkerStyle(1)
        if variable == "mea_mass_reg":
            color = ROOT.RooFit.LineColor(3)
            style_histo = ROOT.RooFit.MarkerStyle(2)

        Roo_histo_1.plotOn(frame, color, name)
        RooBukinPdf.plotOn(frame, color, name)

        cur_mean = Xp.getVal()
        cur_mean_error = np.abs(Xp.getError())

        cur_sig = sP.getVal()
        cur_sig_error =  np.abs(sP.getError())

#        cur_error = cur_sig_error/cur_mean + cur_mean_error*cur_sig/(cur_mean**2)
        cur_error = cur_sig_error/cur_sig
        if variable == "est_mass":
            errors.append(cur_error)

        cur_mean = str(cur_mean)
        cur_sig = str(sP.getVal())
        cur_unc = str(sP.getVal() / Xp.getVal())

        leg.AddEntry(frame.getCurve(variable), variable)
        leg.AddEntry(frame.getCurve(variable), "mean = " + cur_mean)
        leg.AddEntry(frame.getCurve(variable), "sigma = " + cur_sig)
        leg.AddEntry(frame.getCurve(variable), "rel-uncert. = " + cur_unc)

    labeler += 1
    data.Close()
    leg.SetTextSize(0.015)
    c1.cd()
    frame.Draw()
    leg.Draw()
    c1.SaveAs("Test_"+str(labeler)+".pdf")

for value in Chisquare_values:
    value = float(value)
print "Relative uncertainties extracted from fit, [alternating between est_mass & mea_mass_reg]"
print rel_uncertainties
print "And the corresponding errors:"
print errors
plt.plot(Chisquare_values, errors)
plt.show()
a = raw_input()
