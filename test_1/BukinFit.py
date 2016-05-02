import ROOT
import numpy as np
ROOT.gSystem.Load("libRooFit")

bin_number = 35
x_boundary = 40
y_boundary = 180

labeler = 0
some_values = []
for no_jets in ["==2","==3","==4", ">4"]:

    bla = ROOT.RooDataHist()

    labeler += 1
    c1 = ROOT.TCanvas("c1", "c1", 800,800)
    leg = ROOT.TLegend(0.68, 0.65, 0.9, 0.9)
 
    x = ROOT.RooRealVar("x-axis", "x", 75.0 , 165.0)
    Xp = ROOT.RooRealVar("Xp_sgn", "Xp",100.0 , 200.0)
    sP = ROOT.RooRealVar("sP_sgn", "sP", 0.0 , 50.0)
    xi = ROOT.RooRealVar("xi_sgn", "xi", -5.0, 5.0)
    rho1 = ROOT.RooRealVar("rho1_sgn", "rho1", -2.0, 2.0)
    rho2 = ROOT.RooRealVar("rho2_sgn", "rho2", -2.0, 2.0)
    RooBukinPdf = ROOT.RooBukinPdf("buk_pdf_sgn","RooBukinPdf for signal", x, Xp, sP, xi, rho1, rho2)

    frame = x.frame()
    frame.SetName("frame")
    frame.SetTitle("Higgs mass estimate w. a Bukin function fitted to it for nJet " + no_jets)
    
    mean_corrfac_sigma = []
    iteration_counter = 0

    for mass in ["est_mass", "mea_mass", "superposition"]:
        
        data = ROOT.TFile('ChiSquareFits_NeutrinoAdded_Full_only_Flavour5_hJCidx.root')
        Tree = ROOT.gDirectory.Get('ChiSquareFits')
        
        if mass == "mea_mass":
            Tree.Draw(mass+">>h1(" +str(bin_number) + "," + str(x_boundary) + "," + str(y_boundary) + ")", "nJet" + no_jets )
        
        if mass == "est_mass":
            Tree.Draw(mass+">>h1(" +str(bin_number) + "," + str(x_boundary) + "," + str(y_boundary) + ")", "nJet" + no_jets )

        if iteration_counter in [0,1]:        
            h1 = ROOT.gDirectory.Get("h1")
            Roo_histo_1 = ROOT.RooDataHist("higgs_mass_estimate", "estimated higgs mass",  ROOT.RooArgList(x), h1, 1.0)
            RooBukinPdf.fitTo(Roo_histo_1)
            mean_corrfac_sigma.append(Xp.getVal())
            mean_corrfac_sigma.append(sP.getVal())

        if iteration_counter == 2:
            
            #Extracting the correlation factor from the est_mass & mea_mass diagrams
            Tree.Draw("est_mass*(est_mass > 0) + mea_mass*(est_mass == 0)>>est_histo(" +str(bin_number) + "," + str(x_boundary) + "," + str(y_boundary) + ")", "est_mass < 250 && nJet" + no_jets)
            Tree.Draw("mea_mass>>mea_histo(" +str(bin_number) + "," + str(x_boundary) + "," + str(y_boundary) + ")", "est_mass < 250 && nJet" + no_jets)
            est_histo = ROOT.gDirectory.Get("est_histo")
            mea_histo = ROOT.gDirectory.Get("mea_histo")
        
            upperbound_est = mean_corrfac_sigma[0] + mean_corrfac_sigma[1]
            lowerbound_est = mean_corrfac_sigma[0] - mean_corrfac_sigma[1]
            
            upperbound_mea = mean_corrfac_sigma[2] + mean_corrfac_sigma[3]
            lowerbound_mea = mean_corrfac_sigma[2] - mean_corrfac_sigma[3]

            string_2d = "est_mass > " + str(lowerbound_est) + "&& est_mass <" + str(upperbound_est) + "&& mea_mass > " + str(lowerbound_mea) + "&& mea_mass < " + str(upperbound_mea)
            Tree.Draw("est_mass:mea_mass>>two_dim_histo", "est_mass != 0 &&" + string_2d + " && nJet" + no_jets, "COLZ")
            two_dim_histo = ROOT.gDirectory.Get("two_dim_histo")
            corr_fac = two_dim_histo.GetCorrelationFactor()
            mean_corrfac_sigma.append(corr_fac)
            
            #Calculate the factors for the superposition of the est & mea diagrams
            overall_factor = 1.0/(1.0 / mean_corrfac_sigma[1]**2 - 2.0*corr_fac/(mean_corrfac_sigma[1]*mean_corrfac_sigma[3]) + 1.0 / mean_corrfac_sigma[3]**2)
            est_factor = (1.0/mean_corrfac_sigma[1]**2 - corr_fac/(mean_corrfac_sigma[1]*mean_corrfac_sigma[3]))*overall_factor
            mea_factor = (1.0/mean_corrfac_sigma[3]**2 - corr_fac/(mean_corrfac_sigma[1]*mean_corrfac_sigma[3]))*overall_factor
            constant = -(mean_corrfac_sigma[0] - mean_corrfac_sigma[2])*mea_factor*overall_factor

            est_histo.Scale(est_factor)
            mea_histo.Scale(mea_factor)
            est_histo.Add(mea_histo)
            no_bins = int(est_histo.GetNbinsX())

            for i in xrange(no_bins):
                est_histo.AddBinContent(i, constant)

            bla = ROOT.RooDataHist("higgs_mass_estimate", "Weighted superposition of est_mass & mea_mass",  ROOT.RooArgList(x), est_histo, 1.0)
            RooBukinPdf.fitTo(bla)

            
        some_values.append(Xp.getVal())
        some_values.append(sP.getVal())

        color = ROOT.RooFit.LineColor(2)
        name = ROOT.RooFit.Name(mass)

        style_histo = ROOT.RooFit.MarkerStyle(1)

        if mass == "mea_mass":
            color = ROOT.RooFit.LineColor(3)
            style_histo = ROOT.RooFit.MarkerStyle(2)

        if iteration_counter == 2:
            color = ROOT.RooFit.LineColor(4)
            style_histo = ROOT.RooFit.MarkerStyle(4)

        if iteration_counter != 2:
            Roo_histo_1.plotOn(frame, color)
        else:
            bla.plotOn(frame,color)

        cur_mean = str(Xp.getVal())
        cur_sig = str(sP.getVal())

        RooBukinPdf.plotOn(frame, color, name)
        leg.AddEntry(frame.getCurve(mass), mass)
        leg.AddEntry(frame.getCurve(mass), "mean = " + cur_mean)
        leg.AddEntry(frame.getCurve(mass), "sigma = " + cur_sig)
        iteration_counter += 1

    data.Close()
#    ROOT.gDirectory.DeleteAll()
#    ROOT.gDirectory.Delete("c1")
#    ROOT.gDirectory.Delete("frame")
#    ROOT.gDirectory.Delete("higgs_mass_estimate")
#    ROOT.gDirectory.Delete("buk_pdf_sgn")
#    ROOT.gDirectory.Delete("x-axis")

    leg.SetTextSize(0.015)
    c1.cd()
    frame.Draw()
    leg.Draw()
    c1.SaveAs("PNG/OverlayPlot_"+str(labeler)+".pdf")

print "-------------"
print mean_corrfac_sigma    
print "-------------"
print some_values
print "est_mass followed by mea_mass"
