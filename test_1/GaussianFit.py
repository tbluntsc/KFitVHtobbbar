import numpy as np
import ROOT
import matplotlib.pyplot as plt

#Open the raw data (TTree)
RootFile  = ROOT.TFile( 'Egen_to_Ereco_NeutrinosAdded_full_onlyFlavour5.root' )
Tree = ROOT.gDirectory.Get('Egen_to_Ereco')
entries = Tree.GetEntries()

#Define Eta regions
eta_0 = [0.0,1.0]
eta_1 = [1.0,1.5]
eta_2 = [1.5,2.0]
eta_3 = [2.0,2.5]

regions = []
regions.append(eta_0)
regions.append(eta_1)
regions.append(eta_2)
regions.append(eta_3)

HadronFlavours = []
HadronFlavours.append(0)
HadronFlavours.append(5)

#Create ROOT File to write Fits & Histograms to, w a Subdirectory for each Eta region
Saved_Fits = ROOT.TFile( 'Egen_to_Ereco_fits_NeutrinosAdded_full_onlyFlavour5.root', 'RECREATE')
Saved_Fits.cd()

eta_directories = []
no_fits = 36
interval_width = 5
fit_params = np.zeros((no_fits,2,3,len(regions)))

for region in range(len(regions)):
    eta_dir = Saved_Fits.mkdir("eta" + str(region), "Eta region between " + str(regions[region][0]) + " and " + str(regions[region][1]) + ", Interval width = " + str(interval_width))
    eta_directories.append(eta_dir)
    for flavour in HadronFlavours:
        eta_dir.cd()
        flavour_dir = eta_dir.mkdir("Flavour"+str(flavour), "Flavour == "+str(flavour))

for region in range(len(regions)):

    for flavour in HadronFlavours:

        #Go into directory of eta region, then into flavour directory
        eta_directories[region].cd("Flavour"+str(flavour))

        for i in xrange(no_fits):
            #Initialize Interval in which Gaussian Fit will be made
            interval = [(i+4)*interval_width, (i+5)*interval_width]

            #Create the strings that limit the data to the range of the Eta & RecoPt intervals & to the current Hadron Flavour

            string1 = []
            string1.append('RecoPt-GenPt>>h_')
            string1.append(str(i))
            string1 = "".join(string1)

            string2 = []
            string2.append('RecoPt > ')
            string2.append(str(interval[0]))
            string2.append(' && RecoPt < ')
            string2.append(str(interval[1]))
            string2.append(' && abs(RecoEta) > ')
            string2.append(str(regions[region][0]))
            string2.append(' && abs(RecoEta) < ')
            string2.append(str(regions[region][1]))
            string2.append(' && RecoFlavour == ')
            string2.append(str(flavour))
            string2.append(' && GenPt != 0')
            string2 = "".join(string2)
            
            substr = []
            substr.append('h_')
            substr.append(str(i))
            substr = "".join(substr)
    
            #Create Histogram in Histo, fit it to a Gaussian in param and write the Histogram to the root file opened in the start
            test = Tree.Draw(string1, string2)
            Histo = ROOT.gDirectory.Get(substr)

            #Fit gaussian only between Mean +- FWHM of the histogram instead of the whole histogram
            first = Histo.FindFirstBinAbove(Histo.GetMaximum()/2)
            last = Histo.FindLastBinAbove(Histo.GetMaximum()/2)
            fwhm = Histo.GetBinCenter(last) - Histo.GetBinCenter(first)
            mean = Histo.GetMean()

            gaus_range = ROOT.TF1("gaus_range", "gaus", mean-fwhm, mean+fwhm)

            param = Histo.Fit("gaus_range", "S")
            Histo.Write()

        #Go back to the current Eta Directory
        eta_directories[region].cd()

Saved_Fits.ls()
Saved_Fits.Close()
