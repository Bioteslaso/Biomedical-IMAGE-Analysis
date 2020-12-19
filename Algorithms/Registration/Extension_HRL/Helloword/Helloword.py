from __main__ import vtk, qt, ctk, slicer

#
# HelloPython
#

class HelloPython:
  def __init__(self, parent):
    parent.title = "Hello Python"
    parent.categories = ["Examples"]
    parent.dependencies = []
    parent.contributors = ["Jean-Christophe Fillion-Robin (Kitware)",
                           "Steve Pieper (Isomics)",
                           "Sonia Pujol (BWH)"] # replace with "Firstname Lastname (Org)"
    parent.helpText = """
    Example of scripted loadable extension for the HelloPython tutorial.
    """
    parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.,
Steve Pieper, Isomics, Inc., and Sonia Pujol, Brigham and Women's Hospital and was 
partially funded by NIH grant 3P41RR013218-12S1 (NAC) and is part of the National Alliance 
for Medical Image Computing (NA-MIC), funded by the National Institutes of Health through the 
NIH Roadmap for Medical Research, Grant U54 EB005149.""" # replace with organization, grant and thanks.
    self.parent = parent

#
# qHelloPythonWidget
#

class HelloPythonWidget:
  def __init__(self, parent = None):
    if not parent:
      self.parent = slicer.qMRMLWidget()
      self.parent.setLayout(qt.QVBoxLayout())
      self.parent.setMRMLScene(slicer.mrmlScene)
    else:
      self.parent = parent
    self.layout = self.parent.layout()
    if not parent:
      self.setup()
      self.parent.show()

  def setup(self):
    # Instantiate and connect widgets ...

    # Collapsible button
    sampleCollapsibleButton = ctk.ctkCollapsibleButton()
    sampleCollapsibleButton.text = "A collapsible button"
    self.layout.addWidget(sampleCollapsibleButton)

    # Layout within the sample collapsible button
    sampleFormLayout = qt.QFormLayout(sampleCollapsibleButton)

    # HelloWorld button
    # HelloWorld button
    helloWorldButton = qt.QPushButton("Hello world")
    helloWorldButton.toolTip = "Print 'Hello world' in standard ouput."
    sampleFormLayout.addWidget(helloWorldButton)
    helloWorldButton.connect('clicked(bool)', self.onHelloWorldButtonClicked)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Set local var as instance attribute
    self.helloWorldButton = helloWorldButton

  def onHelloWorldButtonClicked(self):
    print "Hello World !"
    qt.QMessageBox.information(slicer.util.mainWindow(),'SlicerPython','HelloWorld!')


