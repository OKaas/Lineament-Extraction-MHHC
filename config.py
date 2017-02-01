# =====================================================================================
# | Original MHHC algorithm for spatial clustering of the line segments
# | as a tool for lineament extraction.
# |
# | Jakub Silhavy, 2016
# | Faculty of Applied Sciences,University of West Bohemia, Pilsen, Czech Republic
# =====================================================================================

workspace = ".\Input\\"
vysledkyDir = workspace + "Vysledky\\"
zdrojDir = workspace + "DEM\\"

shpLinesWS = "shpLines\\"
mergeWS = "merge\\"
temp = "temp\\"
relevantWS = "relevant\\"
bundleWS = "bundle\\"
clearWS =  "clear\\"




#############################
### Relevant Parameters   ###
#############################
splitField = "split"
relevantMergedName = "relevantMerged.shp"
bundleMergedName = "bundleMerged.shp"

parMeaNon = 2
parMeaRel = 4
parMedNon = 2
parMedRel = 4

relevantT = 3

#############################
### Clusters              ###
#############################

memorySaving = 1
optimalStop = 2000
averageMethod = "centroid" # "lw_average" # "aritmetic", "lw_average", "centers", "centroid"
# clusterT = computation threshold - data to delete
clusterT = 4
# countT = display threshold - data to dont visualize - countT renamed to Filter!
filterCount = 4
## plots ##
yMax = 4
radMax = 5

# ziskani cell size of input DEM
def getCellSize(DEM):
  codeSApos = DEM.rfind("_")
  return int(DEM[codeSApos+1:])

# ziskani SA z nazvu rasteru - SA je folder pro vypocet
def getSA(DEM):
  codeSApos = DEM.find("_")
  codeSA =DEM[0:codeSApos]
  if codeSA == "s":
    SA = "Sumava"
  elif codeSA == "lk":
    SA = "Laka"
  elif codeSA == "k":
    SA = "Kozel"
  elif codeSA == "pr":
    SA = "Prasily"
  elif codeSA == "z":
    SA = "Ziar"
  elif codeSA == "z5":
    SA = "Ziar50"
  elif codeSA == "e":
    SA = "Etiopy"
  elif codeSA == "sp":
    SA = "SumavaPaper"
  elif codeSA == "t":
    SA = "Tatry"
  elif codeSA == "zr":
    SA = "ZiarRotace"
  elif codeSA == "d":
    SA = "TestDMR"
  elif codeSA == "z1":
    SA = "ZiarR1"
  elif codeSA == "z3":
    SA = "ZiarR3"
  elif codeSA == "z4":
    SA = "ZiarR4"
  elif codeSA == "c":
    SA = "Canada"
  elif codeSA == "zc":
    SA = "ZiarRc"
  elif codeSA == "o":
    SA = "Dir8"
  elif codeSA == "tr":
    SA = "Trajectory"
  elif codeSA == "hr":
    SA = "Hronska"
  elif codeSA == "a":
    SA = "Artefact"
  else:
    SA = codeSA
  return SA

# ziskani sourceDEM z nazvu rasteru
def getSourceDEM(DEM):
  codeSApos = DEM.find("_")
  codeSourceDEMpos = DEM[codeSApos+1:].find("_") + codeSApos+1
  codeSourceDEM = DEM[codeSApos+1:codeSourceDEMpos]
  if codeSourceDEM == "zm":
    sourceDEM = "ZM50"
  elif codeSourceDEM == "lls":
    sourceDEM = "LLS"
  elif codeSourceDEM == "d":
    sourceDEM = "DMU25"
  elif codeSourceDEM == "as":
    sourceDEM = "ASTER"
  elif codeSourceDEM == "sr":
    sourceDEM = "SRTM"
  else:
    sourceDEM = codeSourceDEM
  return sourceDEM

def getBufferSizeCluster(DEM):
  sourceDEM = getSourceDEM(DEM)
  if sourceDEM == "LLS":
    # for 5 m resolution
    bufferSizeClusterMin = 100 # 30
    bufferSizeClusterMax = 200 # 15
  else:
    # from 30 m resolution up
    bufferSizeClusterMin = 100
    bufferSizeClusterMax = 200
  return [bufferSizeClusterMin, bufferSizeClusterMax]
