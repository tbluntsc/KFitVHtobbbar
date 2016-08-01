import ROOT
import numpy as np
ROOT.gSystem.Load("libRooFit")

labeler = 0
some_values = []
for no_jets in ["==2","==3","==4", ">4"]:

    bla = ROOT.RooDataHist()

    labeler += 1
    c1 = ROOT.TCanvas("c1", "c1", 1000,1500)
    c1.SetRightMargin(0.02)
    c1.SetTopMargin(0.21)
    c1.SetLeftMargin(0.125)
    c1.SetBottomMargin(0.15)

    leg = ROOT.TLegend(0.125, 0.80, 0.98, 0.98)
    leg.SetTextSize(0.9)

    x = ROOT.RooRealVar("x-axis", "x", 70,160)
    Xp = ROOT.RooRealVar("Xp_sgn", "Xp",100.0 , 150.0)
    sP = ROOT.RooRealVar("sP_sgn", "sP", 0.0 , 50.0)
    xi = ROOT.RooRealVar("xi_sgn", "xi", -5.0, 5.0)
    rho1 = ROOT.RooRealVar("rho1_sgn", "rho1", -2.0, 2.0)
    rho2 = ROOT.RooRealVar("rho2_sgn", "rho2", -2.0, 2.0)
    RooBukinPdf = ROOT.RooBukinPdf("buk_pdf_sgn","RooBukinPdf for signal", x, Xp, sP, xi, rho1, rho2)

    frame = x.frame()
    frame.SetName("frame")
    frame.SetTitle("")

    x_axis = frame.GetXaxis()
    x_axis.SetTitle("Mass (GeV)")
    x_axis.SetTitleSize(0.065)
    x_axis.SetTitleOffset(0.88)
    x_axis.SetLabelSize(0.028)

    y_axis = frame.GetYaxis()
    y_axis.SetTitle("Events ( 1 / 2.5 GeV )")
    y_axis.SetTitleSize(0.065)
    y_axis.SetTitleOffset(0.88)
    y_axis.SetLabelSize(0.028)

    mean_corrfac_sigma = []
    iteration_counter = 0

    for mass in ["est_mass", "mea_mass_reg", "mea_mass", "superposition"]:

        data = ROOT.TFile("Current_Root_Files/V21_ChiSquareFitsUpgrade_hJCidx.root")
#        data = ROOT.TFile("Current_Root_Files/Jet_hadronFlavour_independent/V21_ChiSquareFitsUpgrade_hJCidx.root")
        Tree = ROOT.gDirectory.Get('ChiSquareFits')

        if iteration_counter in [0,1,2]:        

            Tree.Draw(mass+">>h1(80,0,200)", "est_mass != 0 && est_mass < 250 && nJet" + no_jets + "&& Chisquare < 5")# && abs(etas) < 1.0")
            h1 = ROOT.gDirectory.Get("h1")
            Roo_histo_1 = ROOT.RooDataHist("higgs_mass_estimate", "estimated higgs mass",  ROOT.RooArgList(x), h1, 1.0)
            RooBukinPdf.fitTo(Roo_histo_1)

            mean_corrfac_sigma.append(Xp.getVal())
            mean_corrfac_sigma.append(sP.getVal())

        if iteration_counter == 3:
            
            #Extracting the correlation factor from the est_mass & mea_mass diagrams
            Tree.Draw("est_mass*(est_mass > 0) + mea_mass_reg*(est_mass == 0)>>est_histo", "est_mass < 250 && nJet" + no_jets)
            Tree.Draw("mea_mass_reg>>mea_histo", "est_mass < 250 && nJet" + no_jets)

            est_histo = ROOT.gDirectory.Get("est_histo")
            mea_histo = ROOT.gDirectory.Get("mea_histo")
        
            upperbound_est = mean_corrfac_sigma[0] + mean_corrfac_sigma[1]
            lowerbound_est = mean_corrfac_sigma[0] - mean_corrfac_sigma[1]
            
            upperbound_mea = mean_corrfac_sigma[2] + mean_corrfac_sigma[3]
            lowerbound_mea = mean_corrfac_sigma[2] - mean_corrfac_sigma[3]

            string_2d = "est_mass > " + str(lowerbound_est) + "&& est_mass <" + str(upperbound_est) + "&& mea_mass_reg > " + str(lowerbound_mea) + "&& mea_mass_reg < " + str(upperbound_mea)
            Tree.Draw("est_mass:mea_mass_reg>>two_dim_histo", "est_mass != 0 &&" + string_2d + " && nJet" + no_jets, "COLZ")
            two_dim_histo = ROOT.gDirectory.Get("two_dim_histo")
            corr_fac = two_dim_histo.GetCorrelationFactor()
            mean_corrfac_sigma.append(corr_fac)
            
            #Calculate the factors for the superposition of the est & mea diagrams
            overall_factor = 1.0/(1.0 / mean_corrfac_sigma[1]**2 - (2.0*corr_fac)/(mean_corrfac_sigma[1]*mean_corrfac_sigma[3]) + 1.0 / mean_corrfac_sigma[3]**2)
            est_factor = (1.0/mean_corrfac_sigma[1]**2 - corr_fac/(mean_corrfac_sigma[1]*mean_corrfac_sigma[3]))*overall_factor
            mea_factor = (1.0/mean_corrfac_sigma[3]**2 - corr_fac/(mean_corrfac_sigma[1]*mean_corrfac_sigma[3]))*overall_factor
            constant = (mean_corrfac_sigma[2] - mean_corrfac_sigma[0])*(-mea_factor)

            Tree.Draw("%f*est_mass + %f*mea_mass_reg + %f>>bla(80,0,200)"%(est_factor, mea_factor, constant), "est_mass != 0 && est_mass < 250 && nJet" + no_jets + "& Chisquare < 5")# && abs(etas) < 1.0")
            test1 = ROOT.gDirectory.Get("bla")

            bla = ROOT.RooDataHist("higgs_mass_estimate", "Weighted superposition of est_mass & mea_mass_reg",  ROOT.RooArgList(x), test1, 1.0)

            RooBukinPdf.fitTo(bla)
            
        some_values.append(Xp.getVal())
        some_values.append(sP.getVal())

        color = ROOT.RooFit.LineColor(2)
        name = ROOT.RooFit.Name(mass)

        style_histo = ROOT.RooFit.MarkerStyle(1)

        if mass == "mea_mass_reg":
            color = ROOT.RooFit.LineColor(3)
            style_histo = ROOT.RooFit.MarkerStyle(2)

        if mass == "mea_mass":
            color = ROOT.RooFit.LineColor(6)
            style_histo = ROOT.RooFit.MarkerStyle(3)

        if iteration_counter == 3:
            color = ROOT.RooFit.LineColor(4)
            style_histo = ROOT.RooFit.MarkerStyle(4)

        if iteration_counter != 3:
            Roo_histo_1.plotOn(frame, color)

        else:
            bla.plotOn(frame,color)
            
        cur_mean = str(round(Xp.getVal(),2))
        cur_sig = str(round(sP.getVal(),2))
        cur_unc = str(round(sP.getVal() / Xp.getVal(),3))

        RooBukinPdf.plotOn(frame, color, name)

        legEntry = leg.AddEntry(frame.getCurve(mass), mass + ": mean = " + cur_mean + ", sigma = " + cur_sig + ", rel_uncert. = " + cur_unc)
        legEntry.SetTextSize(0.0225)

        if iteration_counter ==3:
            legEntry2 = leg.AddEntry(frame.getCurve(mass), "est_factor = " + str(round(est_factor,2)) + ", mea_factor = " + str(round(mea_factor,2)) + ", corr_factor = " + str(round(corr_fac,2)))
            legEntry2.SetTextSize(0.0225)
#        leg.GetListOfPrimitives().First().SetTextSize(0.01)

        iteration_counter += 1

    data.Close()
    c1.cd()
    frame.Draw()
    leg.Draw()
    c1.SaveAs("MakeThisLookBeautiful"+str(labeler)+".pdf")

print "-------------"
print mean_corrfac_sigma    
print "-------------"
print some_values
print "est_mass followed by mea_mass_reg"
