import csv
import glob
import json
import logging
import math
import numpy
import os
import pickle
import random
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *

#
# Dimensionology
#

class Dimensionology(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Dimensionology"
    self.parent.categories = ["Diffusion"]
    self.parent.dependencies = []
    self.parent.contributors = ["Steve Pieper (Isomics Inc.)"]
    self.parent.helpText = """
This module is used to study tracts.
See more information in <a href="https://github.com/SlicerDMRI/SlicerDMRI#Dimensionology">module documentation</a>.
"""
    self.parent.acknowledgementText = """
Developed as part of "HARMONIZING MULTI-SITE DIFFUSION MRI ACQUISITIONS FOR NEUROSCIENTIFIC ANALYSIS ACROSS AGES AND BRAIN DISORDERS" 5R01MH119222.
This file is based on a template originally developed by Jean-Christophe Fillion-Robin, Kitware Inc., Andras Lasso, PerkLab, and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
"""


#
# DimensionologyWidget
#

class DimensionologyWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # Demos Area
    #
    demosCollapsibleButton = ctk.ctkCollapsibleButton()
    demosCollapsibleButton.text = "Parameters"
    self.layout.addWidget(demosCollapsibleButton)

    # Layout within the dummy collapsible button
    self.demosFormLayout = qt.QFormLayout(demosCollapsibleButton)

    self.hspDimensionologyDemoButton = qt.QPushButton("Run HSP SynthSeg Explorer")
    self.demosFormLayout.addWidget(self.hspDimensionologyDemoButton)
    self.hspDimensionologyDemoButton.connect('clicked()', self.hspDimensionologyDemo)

    # no UI right now

    # Add vertical spacer
    self.layout.addStretch(1)

    self.logic = DimensionologyLogic()


  def cleanup(self):
    pass

  def hspDimensionologyDemo(self):
    self.logic.hspDimensionologyDemo()

#
# DimensionologyLogic
#

class DimensionologyLogic(ScriptedLoadableModuleLogic):
  """
  Use parallel axes to explore tract statistics space

  See: http://syntagmatic.github.io/parallel-coordinates/
  https://github.com/BigFatDog/parcoords-es

  Note also experimented with vtk version, but it has fewer features
  and performance is not any better in practice.

  https://vtk.org/Wiki/VTK/Examples/Python/Infovis/ParallelCoordinatesExtraction
  """

  def __init__(self):
    ScriptedLoadableModuleLogic.__init__(self)

    self.directoryPathByID = {}
    self.tractFileNameByLabel = {}
    self.subjects = []
    self.categoricals = []
    self.radiologyReports = {}

  def hspDimensionologyDemo(self):
    """
    Assumes a mount like this:
      sshfs sdp21@door.nmr.mgh.harvard.edu:/space/alta/1/users /mnt/door

      cp /mnt/door/pieper/data/synthsegAll-CSVs.tar.gz /tmp
      (cd /tmp/; tar xfz synthsegAll-CSVs.tar.gz)

      ln -s /mnt/door/pieper/data/screenshotsAll /tmp/Slicer-pieper/

    Launched with something like this:

      ~/Downloads/Slicer-4.13.0-2021-08-13-linux-amd64/Slicer --additional-module-paths ~/slicer4/latest/SlicerDMRI/Modules/Scripted/Dimensionology/
    """

    # rad reports
    print("Getting radiology_reports")
    radPath = "/mnt/door/radiology_reports/Radiology_Table.csv"
    with open(radPath) as radFile:
      csvReader = csv.reader(radFile)
      headers = csvReader.__next__()
      for row in csvReader:
        reportNumber = row[4]
        self.radiologyReports[reportNumber] = row[9]

    # image measurements
    print("Getting image measurements")
    imageRoot = f"/mnt/door/pieper/data/nii-parallel/"
    segRoot = f"/mnt/door/pieper/data/synthsegAll/"
    csvPattern = f"/mnt/door/pieper/data/synthsegAll/*.csv"
    csvPattern = f"/home/pieper/data/synthsegAll-CSVs/*.csv"
    for csvPath in glob.glob(csvPattern):
      csvName = csvPath.split("/")[-1]
      subjectID = csvName[4:-8]
      try:
        with open(csvPath) as csvFile:
          csvReader = csv.reader(csvFile)
          subjectStats = {}
          subjectStats['id'] = subjectID
          headers = csvReader.__next__()
          indices = csvReader.__next__()
          try:
            values = csvReader.__next__()
          except StopIteration:
            continue
          totalVolume = 0.
          for index in range(1,len(values)):
            subjectStats[headers[index]] = values[index]
            if headers[index] != '':
              totalVolume += float(values[index])
          subjectStats['total volume']= totalVolume

          self.subjects.append(subjectStats)
      except FileNotFoundError:
        print(f"Skipping {csvFilePath}")

    dataToPlotString = json.dumps(self.subjects, indent=2)
    categoricalsString = json.dumps(self.categoricals, indent=2)
    radiologyReportsString = json.dumps(self.radiologyReports, indent=2)

    modulePath = os.path.dirname(slicer.modules.dimensionology.path)
    resourceFilePath = os.path.join(modulePath, "Resources", "HSP-ParCoords-template.html")
    html = open(resourceFilePath).read()
    html = html.replace("%%dataToPlot%%", dataToPlotString)
    html = html.replace("%%categoricals%%", categoricalsString)
    html = html.replace("%%radiologyReports%%", radiologyReportsString)

    self.webWidget = slicer.qSlicerWebWidget()
    self.webWidget.size = qt.QSize(1600,1024)
    # self.webWidget.setHtml(html) #; use file so images load
    self.webWidget.show()

    # save for debugging
    htmlPath = slicer.app.temporaryPath+'/data.html'
    open(slicer.app.temporaryPath+'/data.html', 'w').write(html)
    print(f"Saved to {htmlPath}")
    self.webWidget.url = "file://"+htmlPath



  def showHSPBrushedDimension(self, brushedData, minSize=20*1024*1024):
    dataPath = "/mnt/door/pieper/data"
    subjectID = brushedData['brushedSubjectID']
    inputPath = f"{dataPath}/nii-parallel/{subjectID}.nii"
    segPath = f"{dataPath}/synthsegAll/seg-{subjectID}.nii.gz"
    slicer.util.delayDisplay(f"loading {inputPath} and {segPath}", 300)
    if os.path.getsize(inputPath) > minSize:
      slicer.util.loadVolume(inputPath)
      seg = slicer.util.loadSegmentation(segPath)
      #seg.CreateClosedSurfaceRepresentation()
      slicer.util.delayDisplay(f"Displaying {inputPath} and {segPath}", 300)
    else:
      slicer.util.delayDisplay(f"Skipping based on size {inputPath} and {segPath}", 300)

  def saveScreenshot(self, filePath):
    slicer.util.delayDisplay(f"Saving to {filePath}", 300)
    layoutManager = slicer.app.layoutManager()
    threeDWidget = layoutManager.threeDWidget(0)
    pixmap = threeDWidget.parent().grab()
    pixmap.save(filePath)

  def generateSubjectImages(self):
    dataPath = "/mnt/door/pieper/data/screenshotsAll"
    count = 0
    for subject in self.subjects:
      slicer.mrmlScene.Clear(0)
      slicer.util.delayDisplay(f"Processing {subject}", 200)
      try:
        subjectImagePath = f"{dataPath}/{subject['id']}.jpg"
        if os.path.exists(subjectImagePath):
          print(f"Skipping existing {subjectImagePath}")
        elif subject['total volume'] < 100:
          print(f"Skipping due to lack of data")
        else:
          self.showHSPBrushedDimension({'brushedSubjectID': subject['id']})
          self.saveScreenshot(subjectImagePath)
      except:
        print(f"Couldn't save image for {subject['id']}")
      count += 1
      print(f"Completed {count} of {len(self.subjects)}")



#
# DimensionologyTest
#

class DimensionologyTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear()

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_Dimensionology1()

  def test_Dimensionology1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")

    self.delayDisplay('Test passed')

