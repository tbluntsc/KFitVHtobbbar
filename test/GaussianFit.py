import numpy as np
import ROOT
import matplotlib.pyplot as plt

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

#Create ROOT File to write Fits & Histograms to, w a Subdirectory for each Eta region
Saved_Fits = ROOT.TFile( 'Egen_to_Ereco_fits.root', 'RECREATE')
Saved_Fits.cd()
for region in range(len(regions)):
    Saved_Fits.mkdir("eta" + str(region), "Eta region between " + str(regions[region][0]) + " and " + str(regions[region][1]))

#Decomment this if you want to check if the data read from the ROOT file is sensible (it will print every 10000th entry)
#for entry in xrange(entries):
#    ientry = Tree.LoadTree(entry)
#    nb = Tree.GetEntry(entry)
#    if entry%10000 == 0:
#        print "Entry ", entry, "of ", entries
#        print "GenPt value: ", Tree.GenPt
#        print "RecoPt value: ", Tree.RecoPt
#        print "nJet value: ", Tree.nJet

no_fits = 15
interval_width = 10
fit_params = np.zeros((no_fits,2,3,len(regions)))

for region in range(len(regions)):
    #Go into directory of current eta region
    Saved_Fits.cd("eta" + str(region))
    for i in xrange(no_fits):
    #Initialize Interval in which Gaussian Fit will be made
        interval = [(i+2)*interval_width, (i+3)*interval_width]

    #Create the strings that limit the data to the range of the Eta & RecoPt intervals

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
        string2 = "".join(string2)
        print string2

        substr = []
        substr.append('h_')
        substr.append(str(i))
        substr = "".join(substr)
        print substr
    
    #Create Histogram in Histo, fit it to a Gaussian in param and save the three parameters to fit_params (and their errors)
        test = Tree.Draw(string1, string2)
        Histo = ROOT.gDirectory.Get(substr)
        param = Histo.Fit("gaus", "S")
        Histo.Write()

        fit_params[i,0,0,region] = param.Value(0) 
        fit_params[i,0,1,region] = param.Value(1)
        fit_params[i,0,2,region] = param.Value(2)
        fit_params[i,1,0,region] = param.Error(0)
        fit_params[i,1,1,region] = param.Error(1)
        fit_params[i,1,2,region] = param.Error(2)

#print fit_params[:,0,:,0]
#print fit_params[:,0,:,1]

##Create json object that contains the eta regions, number of fits & interval width 
#meta_data = {'eta0': eta_0, 'eta1': eta_1, 'eta2': eta_2, 'eta3': eta_3, 'no_fits': no_fits, 'interval_width': interval_width} 

Saved_Fits.ls()
Saved_Fits.Close()
