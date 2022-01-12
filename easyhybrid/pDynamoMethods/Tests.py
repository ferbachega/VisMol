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
		
		parameters = {"maxIterations":1000,"rmsGradient":1}
		proj.RunSimulation(parameters,"Geometry_Optimization")

		_pattern = "*:LIG.248:C02"

		proj.SphericalPruning(_pattern,25.0)
		proj.SettingFixedAtoms(_pattern,20.0)

		parameters = {"maxIterations":1000,"rmsGradient":0.1}
		proj.RunSimulation(parameters,"Geometry_Optimization")

		proj.PrintSystems()
		proj.SaveProject()
		proj.FinishRun()

	#---------------------------------------------------
	def QCSystemsSetting(self):

		proj= SimulationProject("TIMTest_QCTests")
		proj.LoadSystemFromSavedProject("TIMTest_SetUp.pkl")

		lig = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.248:*")
		glu = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:*")
		his = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:*")

		selections= [ lig, glu, his ]
		SMOmodels = ["am1","am1dphot","pddgpm3","pm3","pm6","rm1"]

		for smo in SMOmodels:
			proj.SetSMOHybridModel( smo, selections, -3, 1 )

		proj.FinishRun()


#============================================================================
if __name__ == "__main__":
	test = Tests()
	#test.SetTIMsytem()
	test.QCSystemsSetting()
