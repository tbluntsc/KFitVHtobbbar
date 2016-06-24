import ROOT
import numpy as np
from matplotlib import pyplot as plt
ROOT.gSystem.Load("libRooFit")

Chisquare_values = np.zeros(20)
for i in xrange(len(Chisquare_values)):
    Chisquare_values[i] = (1 - i*0.05)

print Chisquare_values
a=raw_input()

differences = []
rel_uncertainties_before_mea = []
rel_uncertainties_before_est = []
rel_uncertainties_after_mea = []
rel_uncertainties_after_est = []
errors_mea = []
errors_est = []
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
    frame.SetTitle("Higgs mass estimate w. a Bukin function fitted to it for Chisquare < " +str(Chisquare))

    for variable in ["est_mass", "mea_mass_reg"]:

        data = ROOT.TFile('Current_Root_Files/V21_ChiSquareFitsUpgrade_hJCidx.root')
        Tree = ROOT.gDirectory.Get('ChiSquareFits')

        #Get number of events before restricting Chisquare
        Tree.Draw(variable+">>h1(50)", "est_mass < 250 && est_mass != 0")
        h1 = ROOT.gDirectory.Get("h1")
 
        Roo_histo_before = ROOT.RooDataHist("higgs mass w/o Chisquare restriction", "higgs mass w/o Chisquare restriction", ROOT.RooArgList(x), h1, 1.0)
        counts_before = h1.Integral()
        
        RooBukinPdf.fitTo(Roo_histo_before)        

        cur_mean = str(Xp.getVal())
        cur_sig = str(sP.getVal())
        cur_unc = str(sP.getVal() / Xp.getVal())

        if variable == "est_mass":
            rel_uncertainties_before_est.append(cur_unc)
        else:
            rel_uncertainties_before_mea.append(cur_unc)

        #Create histogram w. restricted Chisquare and get the difference of events 
        Tree.Draw(variable+">>h2(50)", "est_mass < 250 && est_mass != 0  && Chisquare <" +str(Chisquare))
        h2 = ROOT.gDirectory.Get("h2")

        Roo_histo_1 = ROOT.RooDataHist("higgs mass w. Chisquare restriction ", "estimated higgs mass",  ROOT.RooArgList(x), h2, 1.0)
        counts_after = h2.Integral()
        if variable == "est_mass":
            differences.append((counts_before - counts_after)/counts_before)

        RooBukinPdf.fitTo(Roo_histo_1) 
       
        if variable == "est_mass":
            rel_uncertainties_after_est.append(sP.getVal()/Xp.getVal())
        else:
            rel_uncertainties_after_mea.append(sP.getVal()/Xp.getVal())

        cur_mean = Xp.getVal()
        cur_mean_error = np.abs(Xp.getError())

        cur_sig = sP.getVal()
        cur_sig_error =  np.abs(sP.getError())

        cur_error = cur_sig_error/cur_mean + cur_mean_error*cur_sig/(cur_mean**2)

        if variable == "est_mass":
            errors_est.append(cur_error)

        if variable == "mea_mass_reg":
            errors_mea.append(cur_error)

#    labeler += 1
    data.Close()
#    leg.SetTextSize(0.015)
#    c1.cd()
#    frame.Draw()
#    leg.Draw()
#    c1.SaveAs("Test_"+str(labeler)+".pdf")

print ""
print "Relative uncertainties extracted from fit:"
print rel_uncertainties_before_est, "before, est_mass"
print rel_uncertainties_before_mea, "before, mea_mass"
print ""
print "----"
print "Relative uncertainties extracted from fit:"
print rel_uncertainties_after_est, "after, est_mass"
print rel_uncertainties_after_mea, "after, mea_mass"
print ""
print "----"
print "And the corresponding errors:"
print errors_est, "est_mass"
print errors_mea, "mea_mass"
print ""
print "---"
print "And the percentage of lost events:"
print differences

plot_histogram_estErrors = ROOT.TH1F("histo_errors", "Errors for est_mass relative uncertainty plotted vs Chisquare restriction value", len(Chisquare_values), Chisquare_values[len(Chisquare_values)-1],Chisquare_values[0])

plot_histogram_estErrors.SetMarkerSize(3)
plot_histogram_estErrors.SetMarkerStyle(5)

i = 0
for value in Chisquare_values:
    cur_bin = plot_histogram_estErrors.FindBin(value)
    plot_histogram_estErrors.SetBinContent(cur_bin, errors_est[i]*1e+2)

    i += 1

plot_histogram_differences = ROOT.TH1F("histo_diff", "Percentage of counts lost plotted vs Chisquare restriction value", len(Chisquare_values), Chisquare_values[len(Chisquare_values)-1],Chisquare_values[0])

plot_histogram_differences.SetMarkerSize(3)
plot_histogram_differences.SetMarkerStyle(6)

i = 0
for value in Chisquare_values:
    cur_bin = plot_histogram_differences.FindBin(value)
    plot_histogram_differences.SetBinContent(cur_bin, differences[i])
    i += 1

c = ROOT.TCanvas()
c.SetTitle("Errors and Differences")
plot_histogram_differences.SetOption("P")
plot_histogram_differences.Draw()

plot_histogram_estErrors.SetOption("P")
plot_histogram_estErrors.Draw("SAME")

leg = ROOT.TLegend(.55,.5,.9,.7)
leg.AddEntry(plot_histogram_differences, "Differences [10^-2]")
leg.AddEntry(plot_histogram_estErrors, "est Errors")
leg.Draw()

img = ROOT.TImage.Create()
img.FromPad(c)
img.WriteImage("est_errors.png")
