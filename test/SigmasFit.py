import numpy as np
import ROOT
import re

RootFile = ROOT.TFile('Egen_to_Ereco_fits.root')

#Some shenanigans to extract the number of Eta regions and the number of fits in each region from the RootFile
ListOfKeys = RootFile.GetListOfKeys()
first_key = ListOfKeys[0]
name_of_first_key = first_key.GetName()
first_directory = RootFile.GetDirectory(name_of_first_key)
Keys_in_first_Directory = first_directory.GetListOfKeys()

no_fits = Keys_in_first_Directory.GetEntries()
no_regions = ListOfKeys.GetEntries()

#Extract the actual Eta region boundaries from the titles of the directories of the RootFile to save them later
regions = []

for i in xrange(no_regions):
    current_key = ListOfKeys[i]
    title_of_current_key = current_key.GetTitle()
    
    beg_1 = title_of_current_key.find("between") + 8
    end_1 = title_of_current_key.find(" and")
    beg_2 = title_of_current_key.find("and ") + 4

    current_region = []
    current_region.append(np.float32(title_of_current_key[beg_1:end_1]))
    current_region.append(np.float32(title_of_current_key[beg_2:beg_2+3]))
    regions.append(current_region)

#Initialize Dictionary in which the parameters will be saved in the end
result ={}

for i in xrange(no_regions):
    RootFile.cd("eta"+str(i))
    Test = ROOT.TH1F("Sigmas", "Fitted Sigma values for Eta in ["+str(regions[i][0])+","+str(regions[i][1])+"]" , 320000, 15.0, 200.0)

    for fit in xrange(no_fits):
        #Load Histogram h_fit and from it load the Gaussian function fitted to it
        Histo = ROOT.gDirectory.Get("h_"+str(fit))
        myfunc = Histo.GetFunction("gaus")   

        #Parameter 1 and 2 correspond to the position and the width of the Gaussian fit
        current_sigma = myfunc.GetParameter(2)
        current_energy = myfunc.GetParameter(1)

        #Find bin with value of position of the Gaussian and write Gaussian width to it
        bin_number = Test.FindBin(current_energy)
        Test.SetBinContent(bin_number, current_sigma)
    #This Function is not known if defined at the very start of the macro. No idea why.     
    myfunc = ROOT.TF1("myfunc", "[0] + [1]*x + [2]*sqrt(x)", 20, 200)
    params = Test.Fit("myfunc", "S")
    result["eta"+str(i)] = [[params.Value(0), params.Value(1), params.Value(2)],regions[i]]
    
Test.SetMarkerColor(2)
Test.SetMarkerSize(1)
Test.SetMarkerStyle(5)
Test.Draw("P")

print result['eta0']
print result['eta1']
print result['eta2']
print result['eta3']


Bla = raw_input()
#a = Histo.Draw()
RootFile.Close()
