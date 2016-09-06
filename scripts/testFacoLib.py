"""
This script aims at testing each funtion of the facoLib.py file.

qpython 31.08.2016
"""
import facoLib as fl


msg = "testing new function ...\n"


print msg
print fl.combinedRelUncertainty(100, [40., 80], True)

print msg
print fl.combinedUncertainty([40., 80], True)

print msg
x = 4
y = 6
z = 2.1
myDoubleArray = [[0, 1, 2], [x, y, z]]
myHeader = ["this", "stuff", "is great"]
print fl.makeTable("testOutput", myDoubleArray, myHeader, True, "any caption", True)


print msg
print fl.getDictFromJson("WW", "Nu", True)


print msg
print fl.floatToPercent(0.12)
