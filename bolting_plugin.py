"""Abaqus CAE plugin to find and create couplings and wire features connecting
circular edges matching the size and Parts of the two edges selected by the user.

This is intended to make it much more convenient to fasten layers of midsurface
together using connectors when several circular fastener holes are available but
there are no fasteners.

Carl Osterwisch, October 2024
"""

from abaqusGui import *
__version__ = '0.2.2'

class EdgePickingProcedure(AFXProcedure):

    def __init__(self, owner):

        AFXProcedure.__init__(self, owner)

        self.cmd = AFXGuiCommand(self, method='addConnectors', objectName='bolting')
        self.edge1Kw = AFXObjectKeyword(self.cmd, 'edge1', True)
        self.edge2Kw = AFXObjectKeyword(self.cmd, 'edge2', True)

    def getFirstStep(self):
        self.step1 = AFXPickStep(self, self.edge1Kw,
            'Select edge of first bolt hole', AFXPickStep.EDGES,
            numberToPick=ONE)
        return self.step1

    def getNextStep(self, previous):
        if previous == self.step1:
            step = AFXPickStep(self, self.edge2Kw,
                'Select edge of second bolt hole', AFXPickStep.EDGES,
                numberToPick=ONE)
            step.allowRepeatedSelections(False)
            return step

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
        helpUrl='https://github.com/costerwi',
        applicableModules=['Assembly', 'Interaction', 'Load', 'Mesh'],
        description=__doc__,
        )
