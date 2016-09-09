"""
Script that combines various cross section systematic uncertainties into a single one.
These uncertainties are, pdf, scale, and mass assumption.
They are gotten from the following tiwki pages:
1)https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV 
2)https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec
3)https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
4)https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
The output is a sinlge dictionary writen in a json file and that will be used in the UpDownYieldsCalculator.py script.

qpython 06 09 2016

"""

# standard import
import json
import xml.etree.cElementTree as ET


# my import
import facoLib as fl


# dictionary that will containe the down, central and up value of the XS for each sample
CrossSection_dict={}

#
systShifts=["down","central", "up"]


tree = ET.ElementTree(file='../config/ElElV4.xml')

root = tree.getroot()
datasets = root.find('datasets')
print "found  "  + str(len(datasets)) + " datasets"
datasetNames = []
idataset=0


# loop over the dataset with add=1
for d in datasets:
    name = str(d.attrib['name'])
    if d.attrib['add'] == '1' and "Data" not in name and "QCD" not in name and "DY" not in name:
        datasetNames.append(name)
        
        print name
        

        
        # loop over the systShifts
        for systShift in systShifts:
            value = -99


            if systShift == "central" :
                value = 1.0

            # TTjets
            elif "TTJets" in name:
                if systShift == "down" :
                    value = 1 - fl.combinedRelUncertainty(815,[28.61, 34.38, 21.95], False)
                elif systShift == "up":
                    value = 1 + fl.combinedRelUncertainty(815.96, [19.37, 34.38, 22.67], False)
                else:
                    print "error!! ", systShift, "not in ", systShifts
                    
            # tW
            elif "SingleTop_tW" == name:
                if systShift == "down" :
                    value = 1 - fl.combinedRelUncertainty(136.02,[4.57], False)
                elif systShift == "up":
                    value = 1 + fl.combinedRelUncertainty(136.02, [5.40], False)
                else:
                    print "error!! ", systShift, "not in ", systShifts

            # tbarW
            elif "SingleTop_tbarW" == name:
                if systShift == "down" :
                    value = 1 - fl.combinedRelUncertainty(80.95, [3.61], False)
                elif systShift == "up":
                    value = 1 + fl.combinedRelUncertainty(80.95, [4.06], False)
                else:
                    print "error!! ", systShift, "not in ", systShifts

            # Wjets
            elif "WJetsToLNu" == name:
                if systShift == "down" :
                    value = 1 - fl.combinedRelUncertainty(20508.9, [88.2, 770.9], False) # scale, pdf
                elif systShift == "up":
                    value = 1 + fl.combinedRelUncertainty(20508.9, [165.7, 770.9], False)
                else:
                    print "error!! ", systShift, "not in ", systShifts

            # WW
            elif "WW" in name:
                if systShift == "down" :
                    value = 1 - 0.022 
                elif systShift == "up":
                    value = 1 + 0.025
                else:
                    print "error!! ", systShift, "not in ", systShifts

            # WZ
            elif "WZ" in name:
                if systShift == "down" :
                    value = 1 - fl.combinedRelUncertainty(0.106, [0.0036, 0.0050], False) # W+, m(ll) > 40
                elif systShift == "up":
                    value = 1 + fl.combinedRelUncertainty(0.106, [0.0036, 0.0050], False) # the same, sic!
                else:
                    print "error!! ", systShift, "not in ", systShifts


            # ZZ
            elif "ZZ" in name:
                if systShift == "down" :
                    value = 1 - fl.combinedRelUncertainty(0.0349, [0.0011, 0.0016], False)
                elif systShift == "up":
                    value = 1 + fl.combinedRelUncertainty(0.0349, [0.0011, 0.0016], False) # the same, sic!
                else:
                    print "error!! ", systShift, "not in ", systShifts

            # ZG, for now this is factis!!!
            elif "ZG" in name:
                if systShift == "down" :
                    value = 1 - fl.combinedRelUncertainty(0.0349, [0.0011, 0.0016], False)
                elif systShift == "up":
                    value = 1 + fl.combinedRelUncertainty(0.0349, [0.0011, 0.0016], False) # the same, sic!
                else:
                    print "error!! ", systShift, "not in ", systShifts

            # WG, for now this is factis!!!
            elif "WG" in name:
                if systShift == "down" :
                    value = 1 - fl.combinedRelUncertainty(0.0349, [0.0011, 0.0016], False)
                elif systShift == "up":
                    value = 1 + fl.combinedRelUncertainty(0.0349, [0.0011, 0.0016], False) # the same, sic!
                else:
                    print "error!! ", systShift, "not in ", systShifts

            # all the signal
            elif "NP" in name:
                error = -999

                # split in mass, all uncertainties are symetric

                # m=200
                if "m200_" in name:
                    error = fl.combinedRelUncertainty(64.5085, [9.2955], False)
                    if systShift == "down" :
                        value = 1 - error
                    elif systShift == "up":
                        value = 1 + error
                    else:
                        print "error!! ", systShift, "not in ", systShifts

                # m = 300
                if "m300_" in name:
                    error = fl.combinedRelUncertainty(8.51615, [1.18564], False)
                    if systShift == "down" :
                        value = 1 - error
                    elif systShift == "up":
                        value = 1 + error
                    else:
                        print "error!! ", systShift, "not in ", systShifts


                # m = 400
                if "m400_" in name:
                    error = fl.combinedRelUncertainty(1.83537, [0.25142], False)
                    if systShift == "down" :
                        value = 1 - error
                    elif systShift == "up":
                        value = 1 + error
                    else:
                        print "error!! ", systShift, "not in ", systShifts

                # m = 500
                if "m500_" in name:
                    error = fl.combinedRelUncertainty(0.51848, [0.06937], False)
                    if systShift == "down" :
                        value = 1 - error
                    elif systShift == "up":
                        value = 1 + error
                    else:
                        print "error!! ", systShift, "not in ", systShifts

                # m = 600
                if "m600_" in name:
                    error = fl.combinedRelUncertainty(0.174599, [0.023060], False)
                    if systShift == "down" :
                        value = 1 - error
                    elif systShift == "up":
                        value = 1 + error
                    else:
                        print "error!! ", systShift, "not in ", systShifts

                # m = 700
                if "m700_" in name:
                    error = fl.combinedRelUncertainty(0.0670476, [0.0089461], False)
                    if systShift == "down" :
                        value = 1 - error
                    elif systShift == "up":
                        value = 1 + error
                    else:
                        print "error!! ", systShift, "not in ", systShifts

                # m = 800
                if "m800_" in name:
                    error = fl.combinedRelUncertainty(0.0283338, [0.0040152], False)
                    if systShift == "down" :
                        value = 1 - error
                    elif systShift == "up":
                        value = 1 + error
                    else:
                        print "error!! ", systShift, "not in ", systShifts

                # m = 900
                if "m900_" in name:
                    error = fl.combinedRelUncertainty(0.0128895, [0.0019595], False)
                    if systShift == "down" :
                        value = 1 - error
                    elif systShift == "up":
                        value = 1 + error
                    else:
                        print "error!! ", systShift, "not in ", systShifts

                # m = 1000
                if "m1000_" in name:
                    error = fl.combinedRelUncertainty(0.00615134, [0.0010024], False)
                    if systShift == "down" :
                        value = 1 - error
                    elif systShift == "up":
                        value = 1 + error
                    else:
                        print "error!! ", systShift, "not in ", systShifts

                # m = 1100
                if "m1100_" in name:
                    error = fl.combinedRelUncertainty(0.00307413, [0.0005330], False)
                    if systShift == "down" :
                        value = 1 - error
                    elif systShift == "up":
                        value = 1 + error
                    else:
                        print "error!! ", systShift, "not in ", systShifts

                # m = 1200
                if "m1200_" in name:
                    error = fl.combinedRelUncertainty(0.00159844, [0.0002960], False)
                    if systShift == "down" :
                        value = 1 - error
                    elif systShift == "up":
                        value = 1 + error
                    else:
                        print "error!! ", systShift, "not in ", systShifts

                

                    



            # 
            else :
                print "Hum... The name of dataset (" + name + ") does not correspond to any value's formula. You have probably done something wrong and I assume the value is = -99?!"
                print "The value is ", value


                
                

            # write the dict
            CrossSection_dict[name + "_" + systShift]=value
            

            # end of systshifs loop





        

print CrossSection_dict
# save dict in json
with open("jsonFiles/"+"CrossSection"+'.json', 'w') as f:
    json.dump(CrossSection_dict, f)            
                
                        
