import numpy as np
import ROOT
import re

RootFile = ROOT.TFile('Egen_to_Ereco_fits.root')

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
        Test = ROOT.TH1F("Sigmas", "Fitted Sigma values for Eta in ["+str(regions[i][0])+","+str(regions[i][1])+"]" +" and HadronFlavour == "+str(flavour) , 320000, 15.0, 200.0)

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
        myfunc = ROOT.TF1("myfunc", "[0]*x + [1]*sqrt(x) + [2]", 20, 200)
        params = Test.Fit("myfunc", "S")
        result["eta"+str(i),flavour] = [[params.Value(0), params.Value(1), params.Value(2)],regions[i], flavour]

Test.SetMarkerColor(2)
Test.SetMarkerSize(1)
Test.SetMarkerStyle(5)
Test.Draw("P")

print "Hey man check out the eta0 results :"
print result["eta0", 0]
print result["eta0", 5]
print " "
print " "
print "Also eta1 is a special kind of nice :"
print result["eta1", 0]
print result["eta1", 5]
print " "
print " "
print "Have you seen eta2 though? Whew :"
print result["eta2", 0]
print result["eta2", 5]
print " "
print " "
print "At eta3 we have some real gems as well :"
print result["eta3", 0]
print result["eta3", 5]



Bla = raw_input()
#a = Histo.Draw()
RootFile.Close()
