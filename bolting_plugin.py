
from abaqusGui import *

class MyProcedure(AFXProcedure):

    def __init__(self, owner):

        AFXProcedure.__init__(self, owner)

        self.cmd = AFXGuiCommand(self, method='myMethod', objectName='bolting')
        self.edge1Kw = AFXObjectKeyword(self.cmd, 'edge1', True)
        self.edge2Kw = AFXObjectKeyword(self.cmd, 'edge2', True)

    def getFirstStep(self):
        self.step1 = AFXPickStep(self, self.edge1Kw,
            'Select bolt hole edge 1', AFXPickStep.EDGES,
            numberToPick=ONE)
        return self.step1

    def getNextStep(self, previous):
        if previous == self.step1:
            return AFXPickStep(self, self.edge2Kw,
                'Select bolt hole edge 2', AFXPickStep.EDGES,
                numberToPick=ONE)

###########################################################################
# Register the plugins
###########################################################################
toolset = getAFXApp().getAFXMainWindow().getPluginToolset()

toolset.registerGuiMenuButton(
        buttonText='Bolting',
        object=MyProcedure(toolset),
        kernelInitString='import bolting',
        author='Carl Osterwisch',
        #version=__version__,
        #helpUrl=helpUrl,
        applicableModules=['Assembly', 'Interaction', 'Load', 'Mesh'],
        description='',
        )

