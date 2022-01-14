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
		'''
		'''
		proj= SimulationProject("TIMTest_SMOs")
		proj.LoadSystemFromSavedProject("TIMTest_SetUp.pkl")

		lig = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.248:*")
		glu = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:*")
		his = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:*")

		selections= [ lig, glu, his ]
		SMOmodels = ["am1","am1dphot","pddgpm3","pm3","pm6","rm1"]

		#saving qc/mm setup
		proj.SaveProject()
		for smo in SMOmodels:
			proj.SetSMOHybridModel( smo, selections, -3, 1 )

		proj.FinishRun()
	
	#---------------------------------------------------
	def QCDFTBplus(self):
		'''
		'''
		proj= SimulationProject("TIMTest_DFTB")
		proj.LoadSystemFromSavedProject("TIMTest_SetUp.pkl")

		lig = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.248:*")
		glu = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:*")
		his = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:*")

		selections= [ lig, glu, his ]

		proj.SetDFTBsystem(selections, -3, 1 )

		proj.FinishRun()
	
	#---------------------------------------------------
	def QCMMOrca(self):
		'''
		'''
		proj= SimulationProject("TIMTest_ORCA")
		proj.LoadSystemFromSavedProject("TIMTest_SetUp.pkl")

		lig = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.248:*")
		glu = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:*")
		his = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:*")

		selections= [ lig, glu, his ]
		proj.SetOrcaSystem("HF","6-31G*",selections, -3, 1 )
		proj.RunSinglePoint()
		proj.FinishRun()
	
	#---------------------------------------------------
	def QCMMoptimizations(self):
		'''
		'''
		proj= SimulationProject("TIMTest_QCMMopts")
		proj.LoadSystemFromSavedProject("TIMTest_SetUp.pkl")

		lig = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.248:*")
		glu = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:*")
		his = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:*")
		selections= [ lig, glu, his ]
		proj.SetSMOHybridModel( "am1", selections, -3, 1 )
		
		initialCoords = Clone(proj.cSystem.coordinates3)

		#Quasi-Newton muito demorado
		algs = ["ConjugatedGradient",
				"LFBGS",
				"SteepestDescent",
				#"QuasiNewton",
				"FIRE"]
		
		for alg in algs:
			parameters = {"maxIterations":1000,"rmsGradient":0.1,"optmizer":alg}
			proj.RunSimulation(parameters,"Geometry_Optimization")
			proj.cSystem.coordinates3 = initialCoords;

		proj.FinishRun()

	#---------------------------------------------------
	def MD_protocols(self):
		
		proj=SimulationProject("TIMtest_MD")
		proj.LoadSystemFromSavedProject("TIMTest_SetUp.pkl")

		protocols = [
					"heating",
					"equilibration",
					"production"
					]

		for protocol in protocols:
			parameters = {"protocol":protocol,"production_nsteps":2000,"equilibration_nsteps":1000 }
			proj.RunSimulation(parameters,"Molecular_Dynamics")

	#---------------------------------------------------
	def MD_Algs(self):
		'''
		'''
		proj=SimulationProject("TIMtest_MDAlgs")		
		proj.LoadSystemFromSavedProject("TIMTest_SetUp.pkl")

		integrators = [
					"Verlet",
					"LeapFrog",
					"Langevin"
					]

		for integrator in integrators:
			parameters = {"protocol":"production","production_nsteps":2000,"equilibration_nsteps":1000,"MD_method":integrator }
			proj.RunSimulation(parameters,"Molecular_Dynamics")

	#---------------------------------------------------
	def QCMM_MD(self):
		'''
		'''
		proj=SimulationProject("TIMtest_QCMM_MDs")		
		proj.LoadSystemFromSavedProject("TIMTest_SetUp.pkl")

		lig = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.248:*")
		glu = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:*")
		his = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:*")

		selections= [ lig, glu, his ]

		proj.SetSMOHybridModel( "am1", selections, -3, 1 )
		#testing qcmm MD 
		parameters = {"protocol":"production","production_nsteps":2000,"equilibration_nsteps":1000,"MD_method":"LeapFrog" }
		proj.RunSimulation(parameters,"Molecular_Dynamics")

	#---------------------------------------------------
	def QCMM_MDrestricted(self):
		'''
		'''		
		proj=SimulationProject("TIMtest_QCMM_restrictMDs")		
		proj.LoadSystemFromSavedProject("TIMTest_SetUp.pkl")

		lig = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.248:*")
		glu = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:*")
		his = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:*")

		selections= [ lig, glu, his ]

		proj.SetSMOHybridModel( "am1", selections, -3, 1 )
		#testing qcmm MD 

		atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
		atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
		atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")

		'''
		atom4 = AtomSelection.FromAtomPattern(proj.cSystem,"*:OXM.*:O").selection.pop()
		atom5 = AtomSelection.FromAtomPattern(proj.cSystem,"*:OXM.*:H3").selection.pop()
		atom6 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HID.193:NE2").selection.pop()
		'''		
		atomsf = [atom1[0],atom2[0],atom3[0]]

		parameters = {"protocol":"production","production_nsteps":2000,"equilibration_nsteps":1000,"MD_method":"LeapFrog",
					 "atoms":atomsf,"natoms":3, "forceC":100.0,"ndim":1,"MultD1":"true" }
		
		proj.RunSimulation(parameters,"Restricted_Molecular_Dynamics")


	#---------------------------------------------------
	def QCMMScans(self):

		proj=SimulationProject("TIMtest_QCMM_restrictMDs")		
		proj.LoadSystemFromSavedProject("TIMTest_SetUp.pkl")

		lig = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.248:*")
		glu = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:*")
		his = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:*")

		selections= [ lig, glu, his ]

		proj.SetSMOHybridModel( "am1", selections, -3, 1 )
		#testing qcmm MD 

		atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
		atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
		atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")
	
	#---------------------------------------------------
	def UmbrellaSampling(self):
		pass
	
	#---------------------------------------------------
	def ReacCoordSearchers(self):
		pass 
	
	#---------------------------------------------------
	def SMOEnergyRef(self):
		pass
	
	#---------------------------------------------------
	def Thermodynamics(self):
		pass

	#---------------------------------------------------
	def GetOrbitalsInfo(self):
		pass
	
	#---------------------------------------------------
	def Plotting(self):
		pass 
	
	#---------------------------------------------------


#============================================================================
if __name__ == "__main__":
	logFile.Header()
	test = Tests()
	#test.SetTIMsytem()
	#test.QCSystemsSetting()
	#test.QCDFTBplus()
	#test.QCMMOrca()
	#test.QCMMoptimizations()
	#test.MD_protocols()
	#test.MD_Algs()
	#test.QCMM_MD()
	#test.QCMM_MDrestricted()
	test.QCMMScans()
	logFile.Footer()
