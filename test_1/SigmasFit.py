import numpy as np
import ROOT
import re

pow = ROOT.TMath.Power
sqrt = ROOT.TMath.Sqrt

RootFile = ROOT.TFile('Egen_to_Ereco_fits_NeutrinosAdded_full_onlyFlavour5.root')
SigmasFit = ROOT.TFile('Sigmas_Fits_NeutrinoAdded_full_onlyFlavour5.root', 'RECREATE')

#Some shenanigans to extract the number of Eta regions and the number of fits in each region from the RootFile

ListOfKeys = RootFile.GetListOfKeys()
first_key = ListOfKeys[0]
name_of_first_key = first_key.GetName()
first_directory = RootFile.GetDirectory(name_of_first_key)

first_directory.cd()

Keys_in_first_Directory = first_directory.GetListOfKeys()
first_key_in_eta_dir = Keys_in_first_Directory[0]
name_of_first_key_in_eta_dir = first_key_in_eta_dir.GetName()
first_directory_in_eta_dir = first_directory.GetDirectory(name_of_first_key_in_eta_dir)

Keys_in_first_eta_dir = first_directory_in_eta_dir.GetListOfKeys()

no_fits = Keys_in_first_eta_dir.GetEntries()
no_regions = ListOfKeys.GetEntries()

print "no_fits: ", no_fits
print "no_regions: ", no_regions

#Extract the actual Eta region boundaries from the titles of the directories of the RootFile to save them later
regions = []

for i in xrange(no_regions):
    current_key = ListOfKeys[i]
    title_of_current_key = current_key.GetTitle()

    beg_1 = title_of_current_key.find("between") + 8
    end_1 = title_of_current_key.find(" and")
    beg_2 = title_of_current_key.find("and ") + 4
    
    interval_width = title_of_current_key[beg_2+22:beg_2+25]
    interval_width = float(interval_width)
    current_region = []
    current_region.append(np.float32(title_of_current_key[beg_1:end_1]))
    current_region.append(np.float32(title_of_current_key[beg_2:beg_2+3]))
    regions.append(current_region)


#Initialize Dictionary in which the parameters will be saved in the end
result ={}

for i in xrange(no_regions):
    current_eta_dir = RootFile.GetDirectory("eta"+str(i))
    for flavour in [0,5]:
        current_eta_dir.cd("Flavour"+str(flavour))
        Test = ROOT.TH1F("Sigmas_"+str(i)+"_"+str(flavour), "Fitted Sigma values for |Eta| in ["+str(regions[i][0])+","+str(regions[i][1])+"]" +" and HadronFlavour == "+str(flavour) , no_fits, 2.0*interval_width , (no_fits+2)*(interval_width))
        Mean = ROOT.TH1F("Mean_"+str(i)+"_"+str(flavour), "Fitted Mean values for |Eta| in "+str(regions[i][0])+","+str(regions[i][1])+"]" +" and HadronFlavour == "+str(flavour), no_fits, 2.0*interval_width, (no_fits+2)*(interval_width))

        Test.SetMarkerSize(3)
        Test.SetMarkerStyle(5)
        Mean.SetMarkerSize(3)
        Mean.SetMarkerStyle(5)

        for fit in xrange(no_fits):
            #Load Histogram h_fit and from it load the Gaussian function fitted to it
            Histo = ROOT.gDirectory.Get("h_"+str(fit))
            myfunc = Histo.GetFunction("gaus_range")   

            #Parameter 1 and 2 correspond to the position and the width of the Gaussian fit
            current_sigma = myfunc.GetParameter(2)
            current_energy = myfunc.GetParameter(1)

            #Find bin with value of position of the Gaussian and write Gaussian width to it
            pos = ((fit+2)*interval_width + (fit+3)*interval_width)/2.0
            bin_number = Test.FindBin(pos)
            bin_number_mean = Mean.FindBin(pos)
            Test.SetBinContent(bin_number, current_sigma)
            Mean.SetBinContent(bin_number_mean, current_energy)

        #This Function is not known if defined at the very start of the macro. No idea why.     
        sigma_func = ROOT.TF1("sigma_func", "sqrt( [0]*[0]*x*x  + [1]*[1]*x   + [2]*[2])",20,200)
        params = Test.Fit("sigma_func", "S")

        mean_func = ROOT.TF1("mean_func", "[0]*x + [1]", 0, 200)
        params_mean = Mean.Fit("mean_func", "S")

        c = ROOT.TCanvas()
        Test.SetOption("P")
        Test.Draw()
        img = ROOT.TImage.Create()
        img.FromPad(c)
        img.WriteImage("sigma_"+str(i)+"_"+str(flavour)+".png")

        c = ROOT.TCanvas()
        Mean.SetOption("P")
        Mean.Draw()
        img = ROOT.TImage.Create()
        img.FromPad(c)
        img.WriteImage("mean_"+str(i)+"_"+str(flavour)+".png")

        SigmasFit.cd()
        Test.Write()
        Mean.Write()


RootFile.Close()
