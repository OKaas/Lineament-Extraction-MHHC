# =====================================================================================
# | Original MHHC algorithm for spatial clustering of the line segments
# | as a tool for lineament extraction.
# |
# | Jakub Silhavy, 2016
# | Faculty of Applied Sciences, University of West Bohemia, Pilsen, Czech Republic
# =====================================================================================

import autoLin
import config
import tbe
import sys
from multiprocessing import Process, Lock
# setting paths for workspace
workspace = config.workspace
vysledkyDir = config.vysledkyDir
zdrojDir = config.zdrojDir

DEM = sys.argv[1]

result_file_name = 'clustering_result.txt'

result_file = open(result_file_name, 'w')
  
prvni = DEM.find("_r")
druha = DEM.find("_",prvni+1)
rotation = DEM[prvni+2:druha]
rotation = int(rotation)
rotations = [rotation]
DEM = DEM.replace("r%i" % rotation, "dem")
# setup paths for DEM
SA = config.getSA(DEM)
sourceDEM = config.getSourceDEM(DEM)
# time log
timeLogFileName = vysledkyDir + SA + "\\" + sourceDEM + "\\" + "%s_timeLog.txt" % DEM
autoLin.makeDirs(vysledkyDir + SA + "\\" + sourceDEM)
tb = tbe.timeStamp("", "", timeLogFileName)
#tb.log("start")
lock = Lock()
autoLin.optimizedClusterLine(result_file, DEM, rotations, tb, lock)     
#tb.log("finish")
