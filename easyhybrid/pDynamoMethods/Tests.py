#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#FILE = Tests.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#==============================================================================

import os, glob, sys

from commonFunctions import *
import SimulationsPreset 
from LogFile import LogFile

from pBabel                    import *                                     
from pCore                     import *                                     
from pMolecule                 import *                              
from pMolecule.MMModel         import *
from pMolecule.NBModel         import *                                     
from pMolecule.QCModel         import *
from pScientific               import *                                     
from pScientific.Arrays        import *                                     
from pScientific.Geometry3     import *                                     
from pScientific.RandomNumbers import *                                     
from pScientific.Statistics    import *
from pScientific.Symmetry      import *                                     
from pSimulation               import *
from CoreInterface 			   import SimulationProject


ex_path = "/home/igorchem/VisMol/examples/"
timTop  = os.path.join(ex_path,"TIM","7tim.top")
timCrd  = os.path.join(ex_path,"TIM","7tim.crd")

#===========================================================================
class Tests:
	
	#-------------------------------------------------
	def SetTIMsytem(self):
		'''
		Use the methods of the CoreInterface Class to set the Triosephosphate isomerase System
		'''
		proj = SimulationProject("TIMTest_SetUp")
		proj.LoadSystemFromForceField(timTop,timCrd)
		proj.PrintSystems()
		
		co2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.248:C02")
		proj.SphericalPruning(co2,25.0)
		proj.SettingFixedAtoms(co2,20.0)

		parameters = {"optmizer":"LFBGS"}
		proj.RunSimulation(parameters,"Geometry_Optimization")

		proj.SaveProject()
		proj.FinishRun()



	#---------------------------------------------------
	def QCSystemsSetting(self):

		proj= SimulationProject("TIMTest_QCTests")
		proj.LoadSystemFromSavedProject("TIMTest_SetUp.pkl")

		lig = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.248:*")
		glu = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:*")
		his = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:*")

		selections = [ lig, glu, his ]

		proj.SetSMOHybridModel( "am1", selections, -3, 1 )
		proj.FinishRun()


#============================================================================
if __name__ == "__main__":
	test = Tests()
	test.SetTIMsytem()
