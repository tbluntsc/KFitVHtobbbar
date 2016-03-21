import numpy as np
import ROOT
import matplotlib.pyplot as plt

#Open the raw data (TTree)
RootFile  = ROOT.TFile( 'Egen_to_Ereco.root' )
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
Saved_Fits = ROOT.TFile( 'Egen_to_Ereco_fits.root', 'RECREATE')
Saved_Fits.cd()

eta_directories = []

for region in range(len(regions)):
    eta_dir = Saved_Fits.mkdir("eta" + str(region), "Eta region between " + str(regions[region][0]) + " and " + str(regions[region][1]))
    eta_directories.append(eta_dir)
    for flavour in HadronFlavours:
        eta_dir.cd()
        flavour_dir = eta_dir.mkdir("Flavour"+str(flavour), "Flavour == "+str(flavour))
no_fits = 15
interval_width = 10
fit_params = np.zeros((no_fits,2,3,len(regions)))

for region in range(len(regions)):

    for flavour in HadronFlavours:

        #Go into directory of eta region, then into flavour directory
        eta_directories[region].cd("Flavour"+str(flavour))

        for i in xrange(no_fits):
            #Initialize Interval in which Gaussian Fit will be made
            interval = [(i+2)*interval_width, (i+3)*interval_width]

            #Create the strings that limit the data to the range of the Eta & RecoPt intervals & to the current Hadron Flavour

            string1 = []
            string1.append('GenPt>>h_')
            string1.append(str(i))
            string1 = "".join(string1)

            string2 = []
            string2.append('RecoPt > ')
            string2.append(str(interval[0]))
            string2.append(' && RecoPt < ')
            string2.append(str(interval[1]))
            string2.append(' && RecoEta > ')
            string2.append(str(regions[region][0]))
            string2.append(' && RecoEta < ')
            string2.append(str(regions[region][1]))
            string2.append(' && GenFlavour == ')
            string2.append(str(flavour))
            string2 = "".join(string2)

            substr = []
            substr.append('h_')
            substr.append(str(i))
            substr = "".join(substr)
    
            #Create Histogram in Histo, fit it to a Gaussian in param and write the Histogram to the root file opened in the start
            test = Tree.Draw(string1, string2)
            Histo = ROOT.gDirectory.Get(substr)
            param = Histo.Fit("gaus", "S")
            Histo.Write()
        #Go back to the current Eta Directory
        eta_directories[region].cd()

            #fit_params[i,0,0,region] = param.Value(0) 
            #fit_params[i,0,1,region] = param.Value(1)
            #fit_params[i,0,2,region] = param.Value(2)
            #fit_params[i,1,0,region] = param.Error(0)
            #fit_params[i,1,1,region] = param.Error(1)
            #fit_params[i,1,2,region] = param.Error(2)


Saved_Fits.ls()
Saved_Fits.Close()
