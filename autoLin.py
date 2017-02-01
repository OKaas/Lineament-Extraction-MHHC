# =====================================================================================
# | Original MHHC algorithm for spatial clustering of the line segments
# | as a tool for lineament extraction.
# |
# | Jakub Silhavy, 2016
# | Faculty of Applied Sciences, University of West Bohemia, Pilsen, Czech Republic
# =====================================================================================

# tool's library for AutoLin
import shutil
import os
import time
import config
import math
import arcpy
import line_stats
import tbe
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
arcpy.env.overwriteOutput = True
# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")

class Cluster:
  def __init__(self):
    self.myId = []
    self.IDs = []
    self.Ax = []
    self.Ay = []
    self.Bx = []
    self.By = []
    self.A = []
    self.length = []

# class for working with lines - id, length and azimuth
class line:
  def __init__(self, ID, azimuth, length):
    self.ID = ID
    self.azimuth = azimuth
    self.length = length

# setting paths for workspace
workspace = config.workspace
vysledkyDir = config.vysledkyDir
zdrojDir = config.zdrojDir

###########
### GIS ###
###########

def optimizedClusterLine(result_file, DEM, rotations, tb, lock):
  lock.acquire()
  optimalStop = config.optimalStop
  clearWS = config.clearWS
  shpLinesWS = config.shpLinesWS
  bundleWS = config.bundleWS
  mergeWS = config.mergeWS
  demDirs = [bundleWS, clearWS, mergeWS]

  clusterT = config.clusterT

  SA = config.getSA(DEM)
  sourceDEM = config.getSourceDEM(DEM)

  time_count = 0.0

  tb.log("HEre")
  
  # for each rotation in rotations
  for rotation in rotations:

    start_time = time.clock()
    
    rotateDEM = DEM.replace("dem", "r%i" % rotation)
    #tb.log("optimizedClusterLine: %s" % rotateDEM)
    demWS = vysledkyDir + SA + "\\" + sourceDEM + "\\" + rotateDEM + "\\"
    #tb.log("In rotation")
    tb.log(os.path.abspath(demWS + shpLinesWS))
    #tb.log(demWS + shpLinesWS)
    if os.path.exists(os.path.abspath(demWS + shpLinesWS)):
      
      createDEMdirs(demWS, demDirs)

      # first relevant
      inDir = demWS + shpLinesWS
      relevantMergedName = config.relevantMergedName
      relevantT = config.relevantT
      bufferSize = config.getCellSize(DEM)
      relevantMerged = demWS + mergeWS + relevantMergedName

      tb.log(relevantMerged)

      if not os.path.exists(relevantMerged):
        tb.log("relevant")
        relevant(inDir, relevantT, bufferSize, relevantMerged, DEM)

      bundleLines = demWS + mergeWS + "bundleLines.shp"
      tb.log(bundleLines)
      if not os.path.exists(bundleLines):
        result = arcpy.GetCount_management(relevantMerged)
        relevantCount = int(result.getOutput(0))
        relevantCountStart = relevantCount
        stopInterval = relevantCount / 8
        while (relevantCount > optimalStop):
          # cluster line
          relevantMergedB = demWS + clearWS + "RM_%i.shp" % (relevantCount - stopInterval)
          bundleMerged = demWS + bundleWS + "BM_%i.shp" % (relevantCount - stopInterval)
          makeDirs(getDir(relevantMergedB))
          
          tb.log("cluster")
          cluster(relevantMerged, clusterT, stopInterval, relevantMergedB, bundleMerged, DEM)

          # split merged
          tb.log("splitMerged")
          relevantSplit = relevantMergedB
          relevantSplitName = getName(relevantSplit)
          splitDir = getDir(relevantSplit) + "split_%s\\" % relevantSplitName[:-4]
          makeDirs(splitDir)
          splitMerged(relevantSplit, splitDir)

          # relevant
          bufferSize = config.getCellSize(DEM)
          relevantCleared = relevantMergedB[:-4] + "_C.shp"
          #tb.log("relevant")
          relevant(splitDir, relevantT, bufferSize, relevantCleared, DEM)

          # cycle
          relevantMerged = relevantCleared
          result = arcpy.GetCount_management(relevantMerged)
          relevantCount = int(result.getOutput(0))
          # more than 1/2 remains - 1/8 else 1/4
          if (relevantCountStart / 2 + 2) <= relevantCount:
            stopInterval = relevantCount / 8
          else:
            stopInterval = relevantCount / 4
          print relevantCount

        # cluster line
        relevantMergedB = demWS + clearWS + "RM_%i.shp" % (relevantCount - stopInterval)
        bundleMerged = demWS + bundleWS + "BM_%i.shp" % (relevantCount - stopInterval)
        makeDirs(getDir(relevantMergedB))
        
        tb.log("Clustering ...")
        start_time = time.clock()
        cluster(relevantMerged, clusterT, optimalStop, relevantMergedB, bundleMerged, DEM)
        result_file.write("Time: %f \n" % (time.clock() - start_time))
        tb.log("End clustering ...")

        # Create Average Lines of Bundles
        tb.log("createBundleLines")
        createBundleLines(demWS)
  lock.release()

def cluster(relevantMerged, countThreshold, stopInterval, relevantMergedB, bundleMerged, DEM):
  # config:
  [bufferSizeMin, bufferSizeMax] = config.getBufferSizeCluster(DEM)
  azimuthTreshold = config.azimuthTreshold
  # load relevantMerged to in_memory
  relevantMergeLayer = "relevantMerged"
  myMemoryFeature = "in_memory" + "\\" + relevantMergeLayer
  arcpy.CopyFeatures_management(relevantMerged, myMemoryFeature)
  arcpy.MakeFeatureLayer_management(myMemoryFeature, relevantMergeLayer)

  # calculates ID for deleting (FID is recalculated each time after deleting)
  #print("calculates ID for deleting")
  if (arcpy.ListFields(relevantMergeLayer, "ShapeId_1") == []):
    arcpy.AddField_management(relevantMergeLayer, "ShapeId_1", "LONG", 8, 2, "", "", "NULLABLE", "NON_REQUIRED")
  arcpy.CalculateField_management(relevantMergeLayer, "ShapeId_1", "!FID!", "PYTHON_9.3", "#")

  ### for each row in relevantMerged ###
  # cursor for iterate rows in length descend order!
  blueSet = []
  bundleSet = []
  # order cursor from the longest to the shortest line
  rows = arcpy.SearchCursor(relevantMergeLayer, "", "", "", "length D")  # where clause for testing '"ShapeId_1" = 1151'
  # fill the blueset with the rows
  for row in rows:
    blueLine = line(row.ShapeId_1, row.azimuth, row.length)
    blueSet.append(blueLine)
  del rows

  # cleaning
  relevantBackups = []

  ######################################
  # for each line in blueset           #
  ######################################
  for blueLine in blueSet[:stopInterval]:
    myExpression = '"ShapeId_1" =%i' % (blueLine.ID)
    arcpy.SelectLayerByAttribute_management(relevantMergeLayer, "NEW_SELECTION", myExpression)
    noSelected = int(arcpy.GetCount_management(relevantMergeLayer).getOutput(0))
    # if line with ID exists
    if (noSelected == 1):
      # make buffer around blueSHP (bigBuffer for completely within)
      tempBuffer = "in_memory\\tempBuffer"
      # ? Question of buffer size! # dynamic buffer size according to blueLength - 1/10
      blueLength = blueLine.length
      bufferSize = int(blueLength / 10)
      # supremum of buffersize (for extra long paralel lines near together)
      if (bufferSize > bufferSizeMax):
        bufferSize = bufferSizeMax
      # infimum of buffersize (for short lines - small buffer not sufficient)
      if (bufferSize < bufferSizeMin):
        bufferSize = bufferSizeMin
      arcpy.Buffer_analysis(relevantMergeLayer, tempBuffer, "%d Meters" % bufferSize, "FULL", "ROUND", "ALL", "#")
      arcpy.MakeFeatureLayer_management(tempBuffer, "tempBuffer")
      # select all orange in buffer of blue
      # intersect is better but slower - we will see
      arcpy.SelectLayerByLocation_management(relevantMergeLayer, "COMPLETELY_WITHIN", "tempBuffer", "", "NEW_SELECTION")
      noSelected = int(arcpy.GetCount_management(relevantMergeLayer).getOutput(0))
      isBundle = False

      if (noSelected >= countThreshold):
        # create expression +/- azimuthTreshold from blueLine
        blueMin = blueLine.azimuth - azimuthTreshold
        if blueMin < 0:
          blueMin += 180
        blueMax = blueLine.azimuth + azimuthTreshold
        if blueMax > 180:
          blueMax -= 180
        # this condition is useless. Azimuth is always >=0 and <180, after this simplification the myExpression is the same for both cases. The only important thing is to convert extremes to interval <0,180)
        if (blueLine.azimuth < azimuthTreshold) or (blueLine.azimuth > 180 - azimuthTreshold):
          myExpression = '("azimuth" >= %i and "azimuth" < %i) or ("azimuth" > %i and "azimuth" < %i)' % (
          0, blueMax, blueMin, 180)
        else:
          myExpression = '"azimuth" > %i and "azimuth" < %i ' % (blueMin, blueMax)
        ### SELECT THE CLUSTER LINES ###
        arcpy.SelectLayerByAttribute_management(relevantMergeLayer, "SUBSET_SELECTION", myExpression)
        # get count - if < countThreshold do not save, only delete!
        noSelected = int(arcpy.GetCount_management(relevantMergeLayer).getOutput(0))
        if (noSelected >= countThreshold):
          isBundle = True

      if (isBundle):
        # im_memory bundle
        bundle = "in_memory\\line%i" % blueLine.ID
        try:
          arcpy.Buffer_analysis(relevantMergeLayer, bundle, "%d Meters" % 10, "FULL", "ROUND", "ALL", "#")
          # make layer from in_memory bundle
          arcpy.MakeFeatureLayer_management(bundle, "bundle")
          bundleSet.append(bundle)
        except:
          try:
            arcpy.Buffer_analysis(relevantMergeLayer, bundle, "%d Meters" % 12, "FULL", "ROUND", "ALL", "#")
            # make layer from in_memory bundle
            arcpy.MakeFeatureLayer_management(bundle, "bundle")
            bundleSet.append(bundle)
          except:
            continue

        arcpy.AddField_management(bundle, "count", "LONG", 9, "", "", "", "NULLABLE", "NON_REQUIRED")
        arcpy.AddField_management(bundle, "azimuth", "LONG", 9, "", "", "", "NULLABLE", "NON_REQUIRED")
        arcpy.AddField_management(bundle, "length", "LONG", 9, "", "", "", "NULLABLE", "NON_REQUIRED")

        arcpy.CalculateField_management(bundle, "count", noSelected, "PYTHON_9.3", "#")

        lengthList = []
        azimuthList = []

        # compute stats on selection (cluster lines)
        clusterRows = arcpy.SearchCursor(relevantMergeLayer)
        for clusterRow in clusterRows:
          lengthList.append(clusterRow.getValue("length"))
          azimuthList.append(clusterRow.getValue("azimuth"))
        del clusterRows
        # length stats
        [n, mean, std, median, myMin, myMax] = getProperties(lengthList)
        arcpy.CalculateField_management(bundle, "length", "%i" % int(mean + std), "PYTHON_9.3", "#")

        azimuthList.sort()
        azimuthMin = azimuthList[0]
        azimuthMax = azimuthList[n - 1]

        # solve problem with angle numbers!
        # set is on border of azimuths (180-0)
        if ((azimuthMax - azimuthMin) > (2 * azimuthTreshold)):
          # new set - recclassify
          azimuthListPlus = []
          for azimuth in azimuthList:
            if azimuth > (2 * azimuthTreshold):
              azimuthListPlus.append(azimuth - 180)
            else:
              azimuthListPlus.append(azimuth)
          # replace azimuthList
          azimuthList = azimuthListPlus

        # compute azimuth statistics
        [n, mean, std, median, myMin, myMax] = getProperties(azimuthList)
        if mean < 0:
          mean += 180
        arcpy.CalculateField_management(bundle, "azimuth", "%i" % int(mean), "PYTHON_9.3", "#")
      # delete from merged
      arcpy.DeleteFeatures_management(relevantMergeLayer)

      ####################################  E N D   F O R  ###########################################

  #print "backup"
  # a) backup relevantMerged to disk
  arcpy.SelectLayerByAttribute_management(relevantMergeLayer, "CLEAR_SELECTION")
  arcpy.CopyFeatures_management(relevantMergeLayer, relevantMergedB)
  relevantBackups.append(relevantMergedB)
  # write memory bundles to disk
  toMerge = ""
  for bundle in bundleSet:
    toMerge += "in_memory\\%s;" % bundle
  # TODO: don't merge if toMerge is empty !
  try:
    arcpy.Merge_management(toMerge, bundleMerged)
    relevantBackups.append(bundleMerged)
  except Exception, e:
    print toMerge
    print e
  # free memory
  del bundleSet

def relevant(inDir, countThreshold, bufferSize, relevantMerged, DEM):
  tempWS = config.temp
  relevantWS = config.relevantWS
  demDirs = [tempWS, relevantWS]
  createDEMdirs(inDir, demDirs)
  splitField = config.splitField

  inRaster = inDir + tempWS + "outAllR"
  # create sum of rasters
  if not os.path.exists(inRaster):
    getRelevantRaster(inDir, bufferSize, inRaster, DEM)
  #else:
  #  print "Out all raster exist"

  # select only relevant lines according to raster approach (same as in Negatives)
  shp_listDir = os.listdir(inDir)
  for mySHP in shp_listDir:
    if mySHP[-4:] == ".shp":
      inSHP = inDir + mySHP
      lineSHPRelevantsName = mySHP[:-4] + ".shp"
      lineSHPRelevant = inDir + relevantWS + lineSHPRelevantsName
      if not os.path.exists(lineSHPRelevant):
        # zonal with one statistic
        zonalLineRelevant(inSHP, inRaster, bufferSize, countThreshold, lineSHPRelevant, "MEDIAN")

  # # merge to oneSHP
  if not os.path.exists(relevantMerged):
    toMerge = ""
    shp_listDir = os.listdir(inDir + relevantWS)
    for mySHP in shp_listDir:
      if mySHP[-4:] == ".shp":
        inSHP = inDir + relevantWS + mySHP
        # add attribute splitField
        if (arcpy.ListFields(inSHP, splitField) == []):
          arcpy.AddField_management(inSHP, splitField, "TEXT")
        splitName = mySHP[:-4]
        arcpy.CalculateField_management(inSHP, splitField, "'%s'" % splitName, "PYTHON_9.3", "#")
        toMerge += "%s;" % (inSHP)
    arcpy.Merge_management(toMerge, relevantMerged)
  # compute azimuths and lenghts
  try:
    calcStats(relevantMerged)
  except Exception, e:
    print e
  #print "time %s - deleting relevantWS" % (time.strftime("%m%d_%H%M%S"))
  deletePath(inDir + relevantWS)

def createBundleLines(demWS):
  mergeWS = config.mergeWS
  bundleWS = config.bundleWS
  tempWS = config.temp
  demDirs = [tempWS]
  createDEMdirs(demWS, demDirs)

  #print("cluster line algorithm")
  # merge all shp in bundleWS to one to mergeWS - to function in autoLin!
  # merge to oneSHP
  bundleMerged = demWS + mergeWS + config.bundleMergedName
  toMerge = ""
  shp_listDir = os.listdir(demWS + bundleWS)
  for mySHP in shp_listDir:
    if mySHP[-4:] == ".shp":
      toMerge += "%s;" % (demWS + bundleWS + mySHP)
  if len(toMerge) > 0:
    arcpy.Merge_management(toMerge, bundleMerged)
    # createCentroids
    centroidsPoints = demWS + mergeWS + "bundlePoints.shp"
    arcpy.FeatureToPoint_management(bundleMerged, centroidsPoints, "CENTROID")
    # compute average line points
    centroidsLines = demWS + tempWS + "bundleLines.shp"
    centroidsLinesJoin = demWS + mergeWS + "bundleLines.shp"
    newLineList = []
    shapeName = arcpy.Describe(centroidsPoints).shapeFieldName
    rows = arcpy.SearchCursor(centroidsPoints)
    for row in rows:
      arrayLine = arcpy.Array()
      feat = row.getValue(shapeName)
      pnt = feat.firstPoint
      # azimuth is poorly defined in areas around 0 !
      azimuth = row.getValue("azimuth")
      # upper quartile
      length = row.getValue("length")
      dX = math.sin(math.radians(azimuth)) / 2 * length
      dY = math.cos(math.radians(azimuth)) / 2 * length
      startPoint = arcpy.Point(pnt.X - dX, pnt.Y - dY)
      arrayLine.add(startPoint)
      endPoint = arcpy.Point(pnt.X + dX, pnt.Y + dY)
      arrayLine.add(endPoint)
      plyLine = arcpy.Polyline(arrayLine)
      newLineList.append(plyLine)
    del rows
    arcpy.CopyFeatures_management(newLineList, centroidsLines)

    arcpy.env.qualifiedFieldNames = False
    arcpy.MakeFeatureLayer_management(centroidsLines, "centroidsLines")
    arcpy.MakeFeatureLayer_management(centroidsPoints, "centroidsPoints")
    arcpy.AddJoin_management("centroidsLines", "FID", "centroidsPoints", "FID", "KEEP_ALL")
    arcpy.CopyFeatures_management("centroidsLines", centroidsLinesJoin)
  else:
    print "bundleWS is empty"

# input - shp as a single lines, raster with value, buffer size and parameters for selection
# output - directory with separated lines
def zonalLineRelevant(inSHP, inRaster, bufferSize, countThreshold, outSHP, statistic):
  #print "time %s - starting zonal for %s" % (time.strftime("%m%d_%H%M%S"), getName(inSHP))
  tempWS = config.temp
  demDirs = [tempWS]
  inDir = getDir(inSHP)
  createDEMdirs(inDir, demDirs)

  # buffer inSHP
  lineSHPBuffName = getName(inSHP)[:-4] + "_buff%d" % bufferSize
  lineSHPDir = "in_memory\\" + getName(inSHP)[:-4] + "bf_"
  lineSHPBuff = lineSHPDir + lineSHPBuffName
  if not os.path.exists(lineSHPBuff):
    arcpy.Buffer_analysis(inSHP, lineSHPBuff, "%d Meters" % bufferSize, "FULL", "ROUND", "NONE", "#")

  # zonal stats - zones: bufferSHP, inRaster, stat - mean
  if statistic == "MEAN":
    zonalRName = "zMea"
    statField = "fMEAN"
  elif statistic == "MEDIAN":
    zonalRName = "zMed"
    statField = "fMED"

  zonalRMea = lineSHPDir + zonalRName
  if not os.path.exists(zonalRMea):
    arcpy.gp.ZonalStatistics_sa(lineSHPBuff, "FID", inRaster, zonalRMea, statistic, "DATA")
    arcpy.Delete_management(lineSHPBuff)
  # reclassify zone rasters using threshold for flowAcc field (3 classes - positive, unsure, negative)
  # reclassify MEAN and toPolygon
  zonalRMeaRec = lineSHPDir + zonalRName + "%d" % (countThreshold)
  # recclasification condition
  countThreshold -= 0.01
  if not os.path.exists(zonalRMeaRec):
    arcpy.gp.Reclassify_sa(zonalRMea, "VALUE", "0 %d 1;%d 10000000 3" % (countThreshold, countThreshold), zonalRMeaRec,
                           "NODATA")
    arcpy.Delete_management(zonalRMea)
  # toPolygon
  zonalPMeaRec = zonalRMeaRec + "_P"
  if not os.path.exists(zonalPMeaRec):
    arcpy.RasterToPolygon_conversion(zonalRMeaRec, zonalPMeaRec, "SIMPLIFY", "VALUE")
    arcpy.Delete_management(zonalRMeaRec)
  # select and write inSHP, zonalP, outWS
  # cteate new field flowAccMEAN for lineSHP
  if (arcpy.ListFields(inSHP, statField) == []):
    arcpy.AddField_management(inSHP, statField, "LONG", 9, "", "", "", "NULLABLE", "NON_REQUIRED")
  # set null the fields (because 0 is not rewrited!)

  # select and write attributes
  for gridCode in [1, 3]:
    selectAndWrite(inSHP, zonalPMeaRec, gridCode, statField)
  arcpy.Delete_management(zonalPMeaRec)
  # Separating lines - positive and negative
  # A) negative lines > 0
  outDir = getDir(outSHP)
  if not os.path.exists(outDir):
    os.mkdir(outDir)
  # select by attributes: from "lineSHP" using "fTrue"
  myExpression = '%s > 0' % statField
  arcpy.SelectLayerByAttribute_management("lineSHP", "NEW_SELECTION", myExpression)
  arcpy.CopyFeatures_management("lineSHP", outSHP)

def selectAndWrite(lineSHP, zonalP, gridCode, myField):
  # select by attributes: from zonalP using gridcode
  arcpy.MakeFeatureLayer_management (zonalP, "zonalP")
  zonalPFields = arcpy.ListFields("zonalP")
  gridFieldName = "GRIDCODE"
  for zonalField in zonalPFields:
    if "grid" in zonalField.name.lower():
      gridFieldName = zonalField.name
  myExpression = '"%s" = %d' % (gridFieldName, gridCode)
  arcpy.SelectLayerByAttribute_management ("zonalP", "NEW_SELECTION", myExpression)
  # select by location: from lineSHP which have centroid in selected zonalP
  arcpy.MakeFeatureLayer_management(lineSHP, "lineSHP")
  arcpy.SelectLayerByLocation_management("lineSHP", "have_their_center_in", "zonalP", "", "NEW_SELECTION")
  # write to selected lineSHP attribude gridcode to field zonalP
#   myField = "fMED"
  myValue = 1
  if gridCode == 1:
    myValue = -1
  if gridCode == 2:
    myValue = 0
  arcpy.CalculateField_management("lineSHP",myField,myValue,"PYTHON_9.3","#")
  arcpy.SelectLayerByAttribute_management ("zonalP", "CLEAR_SELECTION")
  arcpy.SelectLayerByAttribute_management ("lineSHP", "CLEAR_SELECTION")

def getRelevantRaster(inDir, bufferSize, outAllR, DEM):
  #print "time %s - computing relevant raster" % (time.strftime("%m%d_%H%M%S"))
  zdrojDir = config.zdrojDir
  SA = config.getSA(DEM)
  sourceDEM = config.getSourceDEM(DEM)
  rasterDEM = zdrojDir + SA + "\\" + sourceDEM + "\\" + DEM
  tempWS = config.temp
  demDirs = [tempWS]
  createDEMdirs(inDir, demDirs)
  cellSize = bufferSize
  # list of rasters to merge
  tiffNames = []
  # for each lineaments according to DEM
  shp_listDir = os.listdir(inDir)
  for shpName in shp_listDir:
    if shpName[-4:] == ".shp":
      inSHP = inDir + shpName
      angleStart = shpName.find("_")+1
      angleEnd = shpName.find("_", angleStart)
      rotateAngle = shpName[shpName.rfind("_")+1:-4]
      outLineName = shpName[0:angleEnd] + "_" +rotateAngle
      outLineR = "in_memory\\" + outLineName
      outLineB = "in_memory\\" + outLineName+"B"
      # buffer inSHP
      lineSHPBuffName = shpName[:-4]+"_buff%d" % bufferSize
      lineSHPBuff = "in_memory\\" + lineSHPBuffName
      if not os.path.exists(lineSHPBuff):
        arcpy.Buffer_analysis(inSHP,lineSHPBuff,"%d Meters" % bufferSize,"FULL","ROUND","NONE","#")
      if not os.path.exists(outLineB):
        # create filed, fill with 1 in order to convert vector to raster
        myField = "binary"
        myValue = 1
        # COMPUTING WITH BUFFER instead of line
        inSHP = lineSHPBuff
        arcpy.AddField_management(inSHP, myField, "SHORT", 2, "", "", "", "NULLABLE", "NON_REQUIRED")
        arcpy.CalculateField_management(inSHP,myField,myValue,"PYTHON_9.3","#")
        # polyline to raster, use this field
        if not os.path.exists(outLineR):
          desc = arcpy.Describe(rasterDEM)
          arcpy.env.extent = desc.extent
          arcpy.PolygonToRaster_conversion(inSHP,myField,outLineR,"CELL_CENTER","NONE",cellSize)
          arcpy.env.extent = "default"
        try:
          # reclassify NoData to 0 !
          arcpy.gp.Reclassify_sa(outLineR,"VALUE","1 1;NODATA 0",outLineB,"DATA")
          tiffNames.append(outLineB)
        except Exception, e:
          print e
      arcpy.Delete_management(outLineR)
      arcpy.Delete_management(lineSHPBuff)
  # raster calculator to sum up every raster
  myExpression = ""
  tif_listDir = tiffNames
  for inTIF in tif_listDir:
    myExpression += "\"" + inTIF + "\" + "
  myExpression = myExpression[0:len(myExpression)-2]
  # EXTENT !!!
  arcpy.gp.RasterCalculator_sa(myExpression, outAllR)

  for inTIF in tif_listDir:
    arcpy.Delete_management(inTIF)

def createDEMdirs(demWS, demDirs):
  for demDir in demDirs:
    if not os.path.exists(demWS+demDir):
      os.mkdir(demWS+demDir)

def makeDirs(path):
  if not os.path.exists(path):
    os.makedirs(path)

def getName(inSHP):
   inSHP = inSHP[inSHP.rfind("\\")+1:]
   return inSHP

def getDir(inSHP):
   inSHP = inSHP[0:inSHP.rfind("\\")+1]
   return inSHP

# delete directory with input path
def deletePath(path):
  if os.path.exists(path):
    if os.path.isdir(path):
      #print "time %s - deleting %s" % (time.strftime("%m%d_%H%M%S"), path)
      shutil.rmtree(path)

def splitMerged(relevantMerged, outDir):
  splitField = config.splitField
  # split relevantMerged
  cursor = arcpy.SearchCursor(relevantMerged, "","", splitField)
  # da je do seznamu
  splitNames = []
  for c in cursor:
    splitNames.append(c.getValue(splitField))
  # odstrani duplikaty
  uniqueNames = {}.fromkeys(splitNames).keys()
#   print uniqueNames
  for splitName in uniqueNames:
    outSHP = outDir + splitName+".shp"
#     print "%s = %s" % (splitField, splitName)
    if not os.path.exists(outSHP):
      arcpy.MakeFeatureLayer_management(relevantMerged,splitName, "%s = '%s'" % (splitField, splitName))
      arcpy.CopyFeatures_management(splitName, outSHP)

# lite version of relevant - just merge all shpLines in inDir
def relevantLite(inDir, relevantMerged):
  # merge to oneSHP
  if not os.path.exists(relevantMerged):
    toMerge = ""
    shp_listDir = os.listdir(inDir)
    for mySHP in shp_listDir:
      if mySHP[-4:] == ".shp":
        inSHP = inDir + mySHP
        toMerge+= "%s;" %(inSHP)
    arcpy.Merge_management(toMerge,relevantMerged)
    calcStats(relevantMerged)

# calculates statistics length and azimuth (from north clockwise)
# write fields "azimuth" and "length" to attribute table of inSHP
def calcStats(inSHP):
  if (arcpy.ListFields(inSHP, "azimuth")== []):
    arcpy.AddField_management(inSHP, "azimuth", "FLOAT", 9, 2, "", "", "NULLABLE", "NON_REQUIRED")
  if (arcpy.ListFields(inSHP, "length")== []):
    arcpy.AddField_management(inSHP, "length", "FLOAT", 9, 2, "", "", "NULLABLE", "NON_REQUIRED")

  rows = arcpy.UpdateCursor(inSHP)
  desc = arcpy.Describe(inSHP)
  shapefieldname = desc.ShapeFieldName
  for row in rows:
    feat = row.getValue(shapefieldname)
    pnts = getPoints(feat)
    if (pnts != []):
      A = pnts[0]
      B = pnts[1]
      azimuth = line_stats.lineAzimuth(A,B)
      length = line_stats.lineLength(A,B)
      row.azimuth = azimuth
      row.length = length
      rows.updateRow(row)

# for each input row line returns list of points
def getPoints(feat):
  partnum = 0
  # Step through each part of the feature
  pnts = []
  for part in feat:
      # Step through each vertex in the feature
      for pnt in feat.getPart(partnum):
              pnts.append(pnt)
      partnum += 1
  return pnts

# inSHP - relevantMerged.shp - outSHP - bundleLines_KIV.shp
def clusterKIV(inSHP, resultSHP, method, tb):
  gisExePath = config.gisExePath
  # polarize SHP to TXT
  polarizedTXT = inSHP[:-4]+"_XYAL_P.txt"
  if not os.path.exists(polarizedTXT):
    tb.log("polarize")
    polarize(inSHP, polarizedTXT)
    tb.log("polarize - done")
  # compute KIV cluster algorithm
  resultTXT = inSHP[:-4]+"_result%s.txt" % method # "a" or "b"
  if not os.path.exists(resultTXT):
    X = config.xKIV
    Y = config.yKIV
    A = config.azimuthTreshold
    if method == "a":
      arguments = "%s %i %i %i %s" %(polarizedTXT, X, Y, A, resultTXT)
    elif method == "b":
      # GIS.exe <input file path> <border X> <border Y> <border azimuth> <filter> <output file path>
      arguments = "%s %i %i %i %i %s" %(polarizedTXT, X, Y, A, 0, resultTXT)
    command = "%s %s" % (gisExePath, arguments)
    tb.log("cluster KIV command")
    os.system(command)
    tb.log("cluster KIV command done")

  tb.log("process results")
  myClusters = readClusters(resultTXT)
  clusterTKIV = config.clusterTKIV
  writeAverageSHP(myClusters, resultSHP, clusterTKIV)
  tb.log("process results done")
  # delete temps
  os.remove(polarizedTXT)
  os.remove(resultTXT)

# export coordinates from input SHP to TXT
# polarize lines to directions 0-180 dg
# polarizace nezafunguje na prechodu 179-0 a v jeho okoli - nevhodne treba pro KIV a! a pak nevhodne pro prumerovani
def polarize(inSHP, outTXTName):
  outTXT = open(outTXTName, "w")
  textLog = "ID; A.X; A.Y; B.X; B.Y; A; L\n"
  outTXT.write(textLog)

  rows = arcpy.SearchCursor(inSHP, "","","","length D")

  desc = arcpy.Describe(inSHP)
  shapefieldname = desc.ShapeFieldName
  for row in rows:
    feat = row.getValue(shapefieldname)
    pnts = getPoints(feat)
    A = pnts[0]
    B = pnts[1]

    dY = (A.Y - B.Y)
    dX = (A.X - B.X)
    atan2 = math.atan2(dY,dX)
    alpha = math.degrees(atan2)

    if 90 >= alpha > -90:
      C = A
      A = B
      B = C
    textLog = "%i;%i;%i;%i;%i;%i;%i;\n" % (row.FID, A.X, A.Y, B.X, B.Y, row.azimuth, row.length)
    outTXT.write(textLog)

  outTXT.close()

# add attributes from row to cluster object 1 row = 1 line
def addLine(cluster, row):
  cluster.IDs.append(int(float(row[1])))
  cluster.Ax.append(float(row[2]))
  cluster.Ay.append(float(row[3]))
  cluster.Bx.append(float(row[4]))
  cluster.By.append(float(row[5]))
  cluster.A.append(float(row[6]))
  cluster.length.append(float(row[7]))

# process result file and create list of clusters
def readClusters(resultTXT):
  results = open(resultTXT, "r")
  cluster = Cluster()
  clusterSet = []
  rows = []
  myId = 0

  #   tb.log("read file")
  for result in results:
    rows.append(result)
  results.close()

  #   tb.log("process clusters")
  addLine(cluster, rows[1].replace(",", ".").split(";"))
  cluster.myId = myId
  myId += 1
  for row in rows[2:]:
    row = row.replace(",", ".")
    row = row.split(";")
    if row[0] == "*":
      clusterSet.append(cluster)
      cluster = Cluster()
      cluster.myId = myId
      myId += 1
    addLine(cluster, row)
  clusterSet.append(cluster)
  return clusterSet

def writeAverageSHP(clusterSet, resultSHP, clusterT):
  averageMethod = config.averageMethod
  azimuthTreshold = config.azimuthTreshold
  polylineList = []
  attributeList = []
  countList = []
  for c in clusterSet:
    # filter clusters with insufficient number of lines
    if len(c.IDs) >= clusterT:
      # if cluster is on the border 179-0-1 (suppose that azimuth threshold has been applied)
      if (max(c.A) - min(c.A)) > 2 * azimuthTreshold:
        # for every lines with A < 0 - switch start to end
        for i in range(0, len(c.A), 1):
          if c.A[i] > 2 * azimuthTreshold:
            Cx = c.Ax[i]
            Cy = c.Ay[i]
            c.Ax[i] = c.Bx[i]
            c.Ay[i] = c.By[i]
            c.Bx[i] = Cx
            c.By[i] = Cy

      if (averageMethod == "aritmetic"):
        # aritmetic average of coordinates
        Ax = (sum(c.Ax) / len(c.Ax))
        Ay = (sum(c.Ay) / len(c.Ay))
        Bx = (sum(c.Bx) / len(c.Bx))
        By = (sum(c.By) / len(c.By))
      elif (averageMethod == "centers"):
        # preserve only centers as cluster's representants
        Ax = c.Ax[0]
        Ay = c.Ay[0]
        Bx = c.Bx[0]
        By = c.By[0]
      elif (averageMethod == "lw_average"):
        sumAx = 0
        sumAy = 0
        sumBx = 0
        sumBy = 0
        for i in range(0, len(c.Ax), 1):
          length = c.length[i]
          sumAx += c.Ax[i] * length
          sumAy += c.Ay[i] * length
          sumBx += c.Bx[i] * length
          sumBy += c.By[i] * length
        sumLength = sum(c.length)
        Ax = sumAx / sumLength
        Ay = sumAy / sumLength
        Bx = sumBx / sumLength
        By = sumBy / sumLength
      elif (averageMethod == "centroid"):
        sumAx = 0
        sumAy = 0
        sumBx = 0
        sumBy = 0
        for i in range(0, len(c.Ax), 1):
          length = c.length[i]
          sumAx += c.Ax[i] * length
          sumAy += c.Ay[i] * length
          sumBx += c.Bx[i] * length
          sumBy += c.By[i] * length
        sumLength = sum(c.length)
        Ax = sumAx / sumLength
        Ay = sumAy / sumLength
        Bx = sumBx / sumLength
        By = sumBy / sumLength
        cX = (Ax + Bx) / 2
        cY = (Ay + By) / 2
      arrayLine = arcpy.Array()
      if (averageMethod == "centroid"):
        azimuthList = c.A
        azimuthList.sort()
        azimuthMin = azimuthList[0]
        azimuthMax = azimuthList[-1]
        # solve problem with angle numbers!
        # set is on border of azimuths (180-0)
        if ((azimuthMax - azimuthMin) > (2 * azimuthTreshold)):
          # new set - recclassify
          azimuthListPlus = []
          for azimuth in azimuthList:
            if azimuth > (2 * azimuthTreshold):
              azimuthListPlus.append(azimuth - 180)
            else:
              azimuthListPlus.append(azimuth)
          # replace azimuthList
          azimuthList = azimuthListPlus

        # compute azimuth statistics
        azimuthStats = getProperties(azimuthList)
        azimuth = azimuthStats[1]
        if azimuth < 0:
          azimuth += 180
        lengthStats = getProperties(c.length)
        # average length
        length = lengthStats[1]
        dX = math.sin(math.radians(azimuth)) / 2 * length
        dY = math.cos(math.radians(azimuth)) / 2 * length
        startPoint = arcpy.Point(cX - dX, cY - dY)
        arrayLine.add(startPoint)
        endPoint = arcpy.Point(cX + dX, cY + dY)
        arrayLine.add(endPoint)
      else:
        startPoint = arcpy.Point(Ax, Ay)
        arrayLine.add(startPoint)
        endPoint = arcpy.Point(Bx, By)
        arrayLine.add(endPoint)
      plyLine = arcpy.Polyline(arrayLine)
      polylineList.append(plyLine)
      attributeList.append(c.IDs[0])
      countList.append(len(c.IDs))
  if not polylineList == []:
    arcpy.CopyFeatures_management(polylineList, resultSHP)
    countField = "count"
    if (arcpy.ListFields(resultSHP, countField) == []):
      arcpy.AddField_management(resultSHP, countField, "LONG", 8, 2, "", "", "NULLABLE", "NON_REQUIRED")

    # update cursor - fill attributes
    rows = arcpy.UpdateCursor(resultSHP)
    i = 0
    for row in rows:
      row.Id = attributeList[i]
      row.count = countList[i]
      rows.updateRow(row)
      i += 1
    del rows
  else:
    print "No clusters created!"

# get statistic properties of input set
def getProperties(mySet):
  mySet.sort()
  n = len(mySet)
  myMin = mySet[0]
  myMax = mySet[n - 1]
  mean = sum(mySet) / len(mySet)
  std = math.sqrt(sum((x - mean) ** 2 for x in mySet) / n)
  medianIndex = int(0.5 * (n + 1))
  median = mySet[medianIndex]
  #   print " n: %i \n mean: %0.2f \n std: %0.2f \n median: %0.2f \n min: %0.2f \n max: %0.2f \n " % (n, mean, std, median, myMin, myMax)
  return [n, mean, std, median, myMin, myMax]
