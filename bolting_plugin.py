"""Generate coupling connections between selected hole edges,
find and repeat for all similar edges
"""

from abaqusGui import *
__version__ = '0.1.0'

class EdgePickingProcedure(AFXProcedure):

    def __init__(self, owner):

        AFXProcedure.__init__(self, owner)

        self.cmd = AFXGuiCommand(self, method='myMethod', objectName='bolting')
        self.edge1Kw = AFXObjectKeyword(self.cmd, 'edge1', True)
        self.edge2Kw = AFXObjectKeyword(self.cmd, 'edge2', True)

    def getFirstStep(self):
        self.step1 = AFXPickStep(self, self.edge1Kw,
            'Select edge of first bolt hole', AFXPickStep.EDGES,
            numberToPick=ONE)
        return self.step1

    def getNextStep(self, previous):
        if previous == self.step1:
            return AFXPickStep(self, self.edge2Kw,
                'Select edge of second bolt hole', AFXPickStep.EDGES,
                numberToPick=ONE)

    def getLoopStep(self):
        return self.step1  # loop until stopped

###########################################################################
# Register the plugins
###########################################################################
toolset = getAFXApp().getAFXMainWindow().getPluginToolset()

toolset.registerGuiMenuButton(
        buttonText='Bolting',
        object=EdgePickingProcedure(toolset),
        kernelInitString='import bolting',
        author='Carl Osterwisch',
        version=__version__,
        #helpUrl=helpUrl,
        applicableModules=['Assembly', 'Interaction', 'Load', 'Mesh'],
        description=__doc__,
        )
