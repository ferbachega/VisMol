#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = Tests.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#==============================================================================

import os, glob, sys
#------------------------------------------------------
from commonFunctions import *
import SimulationsPreset 
#------------------------------------------------------
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
from PotentialOfMeanForce import *
#-------------------------------------------------------
ex_path = "/home/igorchem/VisMol/examples/"
timTop  = os.path.join(ex_path,"TIM","7tim.top")
timCrd  = os.path.join(ex_path,"TIM","7tim.crd")

#*************************************************************************
def SetTIMsytem(self):
	'''
	Use the methods of the CoreInterface Class to set the Triosephosphate isomerase System
	'''
	proj = SimulationProject("TIMTest_SetUp")
	proj.LoadSystemFromForceField(timTop,timCrd)
	
	parameters_a = {"maxIterations":1000,"rmsGradient":1}
	proj.RunSimulation(parameters_a,"Geometry_Optimization")

	_pattern = "*:LIG.248:C02"

	proj.SphericalPruning(_pattern,25.0)
	proj.SettingFixedAtoms(_pattern,20.0)

	parameters_b = {"maxIterations":1000,"rmsGradient":0.1}
	proj.RunSimulation(parameters_b,"Geometry_Optimization")

	proj.PrintSystems()
	proj.SaveProject()
	proj.FinishRun()		

#=================================================================
def QCSystemsSetting(self):
	'''
	'''
	#===========================================
	#TIM
	
	proj= SimulationProject("TIMTest_SMOs")
	proj.LoadSystemFromSavedProject("TIMTest_SetUp.pkl")

	lig = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.248:*")
	glu = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:*")
	his = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:*")

	selections= [ lig, glu, his ]
	SMOmodels = ["am1","am1dphot","pddgpm3","pm3","pm6","rm1"]

	#saving qc/mm setup
	for smo in SMOmodels:
		proj.SetSMOHybridModel( smo, selections, -3, 1 )

	proj.SaveProject()
	proj.FinishRun()
		
#===================================================================
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

#===================================================================
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

#=================================================================
def QCMMoptimizations(self):
	'''
	'''
	#----------------------------------------------
	# TIM
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
			"LFBGS"             ,
			"SteepestDescent"   ,
			"FIRE"              ]
	
	#problems saving trajectory! But geometry opt working
	for alg in algs:
		parameters = {
						"maxIterations":1000,
						"rmsGradient":0.1,
						"optmizer":alg,
						#"save_pdb":"true",
						#"save_traj":"true"
					}
		proj.RunSimulation(parameters,"Geometry_Optimization")
		proj.cSystem.coordinates3 = initialCoords;

	proj.SaveProject()
	proj.FinishRun()
	#*************************************************************
	
#===================================================================
def MD_protocols(self):
	'''
	'''
	proj=SimulationProject("TIMtest_MD")
	proj.LoadSystemFromSavedProject("TIMTest_SetUp.pkl")

	protocols = [
				"heating",
				"equilibration",
				"production"
				]

	for protocol in protocols:
		parameters = {"protocol":protocol,"production_nsteps":2000,"equilibration_nsteps":1000,"MD_method":"Verlet" }
		proj.RunSimulation(parameters,"Molecular_Dynamics")

	proj.FinishRun()

#======================================================================
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

	proj.FinishRun()

#=======================================================================
def QCMM_MD(self):
	'''
	'''
	proj=SimulationProject("TIMtest_QCMM_MDs")		
	proj.LoadSystemFromSavedProject("TIMTest_QCMMopts.pkl")
	#testing qcmm MD 
	parameters = {"protocol":"production","production_nsteps":10000,"equilibration_nsteps":5000,"MD_method":"LeapFrog" }
	proj.RunSimulation(parameters,"Molecular_Dynamics")
	proj.FinishRun()

#========================================================================
def QCMM_MDrestricted(self):
	'''
	'''		
	proj=SimulationProject("TIMtest_MM_restrictMDs")		
	proj.LoadSystemFromSavedProject("TIMTest_SetUp.pkl")		
	#testing qcmm MD 
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")

	atomsf = [ atom1[0], atom2[0], atom3[0] ]

	atom6 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:O06")
	atom5 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")
	#1D Simple distance
	atomss = [ atom4[0], atom5[0], atom6[0] ]

	parameters = {"protocol":"production","production_nsteps":50000,"equilibration_nsteps":10000,"MD_method":"LeapFrog",
				 "atoms_M1":atomsf,"atoms_M2":atomss, "forceC":100.0,"ndim":2,"MC_RC1":True,"MC_RC2":True }
	
	proj.RunSimulation(parameters,"Restricted_Molecular_Dynamics")
	proj.FinishRun()
	#------------------------------------------------------------------
	proj=SimulationProject("TIMtest_QCMM_restrictMDs")		
	proj.LoadSystemFromSavedProject("TIMTest_QCMMopts.pkl")		
	#testing qcmm MD 
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")

	atomsf = [ atom1[0], atom2[0], atom3[0] ]

	atom6 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:O06")
	atom5 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")
	#1D Simple distance
	atomss = [ atom4[0], atom5[0], atom6[0] ]

	parameters = {"protocol":"production","production_nsteps":10000,"equilibration_nsteps":5000,"MD_method":"LeapFrog",
				 "atoms_M1":atomsf,"atoms_M2":atomss, "forceC":100.0,"ndim":2,"MC_RC1":True,"MC_RC2":True }
	
	proj.RunSimulation(parameters,"Restricted_Molecular_Dynamics")
	proj.FinishRun()

#========================================================================
def QCMMScans(self):
	'''
	'''
	proj=SimulationProject("TIMtest_QCMM_Scans1DsimpleDistance")		
	proj.LoadSystemFromSavedProject("TIMTest_QCMMopts.pkl")

	#testing qcmm MD 

	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")
	
	#1D Simple distance
	atomsf = [ atom1[0], atom2[0] ] 

	parameters = { 
					'ATOMS_RC1':atomsf,
					'dincre_RC1':0.1,
					"nSteps_RC1":16,
					"ndim":1,
					"MC_RC1":"true"
				}

	proj.RunSimulation(parameters,"Relaxed_Surface_Scan")		
	proj.FinishRun()

	#===============================================================
	projB=SimulationProject("TIMtest_QCMM_Scans1DmultipleDistance")		
	projB.LoadSystemFromSavedProject("TIMTest_QCMMopts.pkl")
	
	atomsf = [ atom1[0], atom2[0], atom3[0] ] 

	parameters = { 
					'ATOMS_RC1':atomsf,
					'dincre_RC1':0.1,
					"nSteps_RC1":16,
					"ndim":1,
					"MC_RC1":"true"
				}

	projB.RunSimulation(parameters,"Relaxed_Surface_Scan")		
	projB.FinishRun()
	
#=======================================================================
def QCMMScans2D(self):
	
	
	proj=SimulationProject("TIMtest_QCMM_Scans2DsimpleDistance")		
	proj.LoadSystemFromSavedProject("TIMTest_QCMMopts.pkl")
	
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")
	
	atom6 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:O06")
	atom5 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")

	#1D Simple distance
	atomsf = [ atom1[0],atom2[0] ] 
	atomss = [ atom4[0], atom5[0] ]

	parameters = {"maxIterations":1000,"rmsGradient":0.1}
	proj.RunSimulation(parameters,"Geometry_Optimization")

	parameters = { 'ATOMS_RC1':atomsf	  ,
				   'ATOMS_RC2':atomss	  ,
				   'dincre_RC1':0.1 	  ,
				   'dincre_RC2':0.1       , 
				   "nSteps_RC1":6		  ,
				   "nSteps_RC2":6 		  , 
				   "ndim": 2 			  ,
				   "NmaxThreads":        8,
				 }

	proj.RunSimulation(parameters,"Relaxed_Surface_Scan")		
	proj.FinishRun()

	#-----------------------------------------------------------
	projb=SimulationProject("TIMtest_QCMM_Scans2DmixedDistance")		
	projb.LoadSystemFromSavedProject("TIMTest_QCMMopts.pkl")

	atomsf = [ atom1[0], atom2[0], atom3[0] ] 
	atomss = [ atom4[0], atom5[0] ]

	parameters = { 'ATOMS_RC1':atomsf	  ,
				   'ATOMS_RC2':atomss	  ,
				   'dincre_RC1':0.1 	  ,
				   'dincre_RC2':0.1       , 
				   "nSteps_RC1":6		  ,
				   "nSteps_RC2":6 		  , 
				   "ndim": 2 			  ,
				   "MC_RC1":		"true",
				   "NmaxThreads":        8,
				 }

	projb.RunSimulation(parameters,"Relaxed_Surface_Scan")		
	projb.FinishRun()

	#-----------------------------------------------------------
	
	projc=SimulationProject("TIMtest_QCMM_Scans2DmultipleDistance")		
	projc.LoadSystemFromSavedProject("TIMTest_QCMMopts.pkl")

	atomsf = [ atom1[0], atom2[0], atom3[0] ] 
	atomss = [ atom4[0], atom5[0], atom6[0] ]

	parameters = {"maxIterations":1000,"rmsGradient":0.1}
	projc.RunSimulation(parameters,"Geometry_Optimization")

	parameters = { 'ATOMS_RC1':atomsf	  ,
				   'ATOMS_RC2':atomss	  ,
				   'dincre_RC1':0.1 	  ,
				   'dincre_RC2':0.1       , 
				   "nSteps_RC1":6		  ,
				   "nSteps_RC2":6 		  , 
				   "ndim": 2 			  ,
				   "MC_RC1":		"true",
				   "MC_RC2":		"true",
				   "NmaxThreads":        8,
				 }

	projc.RunSimulation(parameters,"Relaxed_Surface_Scan")		
	projc.FinishRun()
	
#===================================================================
def UmbrellaSampling1D(self):
	''''
	'''
	#--------------------------------------------------------------
	proj=SimulationProject("TIMtest_US1Dsimpledistance")		
	proj.LoadSystemFromSavedProject("TIMTest_QCMMopts.pkl")

	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")
	
	atomsf = [ atom1[0], atom2[0] ] 
	
	_path = os.path.join(os.getcwd(),"TIMtest_QCMM_Scans1DsimpleDistance","ScanTraj.ptGeo")

	parameters = { 'ATOMS_RC1':atomsf			,
				   "ndim": 1 					,
				   "samplingFactor":2000 		,
				   "equilibration_nsteps":2000	,
				   "production_nsteps":5000		,
				   "trjFolder":_path 			,
				   "MD_method":"LeapFrog"		,
				   "MC_RC1":"true"				,
				   "NmaxThreads":8 				
				 }

	proj.RunSimulation(parameters,"Umbrella_Sampling")
	proj.FinishRun()
	
	#-----------------------------------------------------------
	projB=SimulationProject("TIMtest_US1DmultipleDistance")		
	projB.LoadSystemFromSavedProject("TIMTest_QCMMopts.pkl")
	
	
	atomsf = [ atom1[0], atom2[0], atom3[0] ] 
	_path = os.path.join(os.getcwd(),"TIMtest_QCMM_Scans1DmultipleDistance","ScanTraj.ptGeo")

	parameters = { 'ATOMS_RC1':atomsf			,
				   "ndim": 1 					,
				   "samplingFactor":2000 		,
				   "equilibration_nsteps":2000	,
				   "production_nsteps":5000		,
				   "trjFolder":_path 			,
				   "MD_method":"LeapFrog"		,
				   "MC_RC1":"true"				,
				   "NmaxThreads":8 				
				 }

	projB.RunSimulation(parameters,"Umbrella_Sampling")
	projB.FinishRun()


#=======================================================================
def UmbrellaSampling2D(self):
	'''
	'''
	proj=SimulationProject("TIMtest_US2DsimpleDistance")		
	proj.LoadSystemFromSavedProject("TIMTest_QCMMopts.pkl")

	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")
	
	atom6 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:O06")
	atom5 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")

	atomsf = [ atom1[0], atom2[0] ] 
	atomss = [ atom4[0], atom5[0] ]

	_path = os.path.join(os.getcwd(),"TIMtest_QCMM_Scans2DsimpleDistance","ScanTraj.ptGeo")

	parameters = { 'ATOMS_RC1':atomsf			,
				   "ATOMS_RC2":atomss			,
				   "ndim": 2 					,
				   "samplingFactor":200 		,
				   "equilibration_nsteps":1000	,
				   "production_nsteps":2000		,
				   "trjFolder":_path 			,
				   "MD_method":"LeapFrog"		,
				   "MC_RC1":"true"				,
				   "MC_RC2":"true"				,
				   "NmaxThreads":8 				
				 }
	
	proj.RunSimulation(parameters,"Umbrella_Sampling")
	proj.FinishRun()

	#----------------------------------------------------------
	projb=SimulationProject("TIMtest_Umbrella2DmixedDistance")		
	projb.LoadSystemFromSavedProject("TIMTest_QCMMopts.pkl")

	atomsf = [ atom1[0], atom2[0], atom3[0] ] 
	atomss = [ atom4[0], atom5[0] ]

	_path = os.path.join(os.getcwd(),"TIMtest_QCMM_Scans2DmixedDistance","ScanTraj.ptGeo")

	parameters = { 'ATOMS_RC1':atomsf			,
				   "ATOMS_RC2":atomss			,
				   "ndim": 2 					,
				   "samplingFactor":200 		,
				   "equilibration_nsteps":1000	,
				   "production_nsteps":2000		,
				   "trjFolder":_path 			,
				   "MD_method":"LeapFrog"		,
				   "MC_RC1":"true"				,
				   "MC_RC2":"true"				,
				   "NmaxThreads":8 				
				 }
	
	projb.RunSimulation(parameters,"Umbrella_Sampling")
	projb.FinishRun()
	#----------------------------------------------------------
	projc=SimulationProject("TIMtest_Umbrella2DmultipleDistance")		
	projc.LoadSystemFromSavedProject("TIMTest_QCMMopts.pkl")

	atomsf = [ atom1[0], atom2[0], atom3[0] ] 
	atomss = [ atom4[0], atom5[0], atom6[0] ]

	_path = os.path.join(os.getcwd(),"TIMtest_QCMM_Scans2DmultipledDistance","ScanTraj.ptGeo")

	parameters = { 'ATOMS_RC1':atomsf			,
				   "ATOMS_RC2":atomss			,
				   "ndim": 2 					,
				   "samplingFactor":200 		,
				   "equilibration_nsteps":1000	,
				   "production_nsteps":2000		,
				   "trjFolder":_path 			,
				   "MD_method":"LeapFrog"		,
				   "MC_RC1":"true"				,
				   "MC_RC2":"true"				,
				   "NmaxThreads":8 				
				 }
	
	projc.RunSimulation(parameters,"Umbrella_Sampling")
	projc.FinishRun()

#=====================================================
def PotentialOfMeanField1D(self):
	'''
	'''
	proj=SimulationProject("TIMtest_PMF1D")		
	proj.LoadSystemFromSavedProject("TIMTest_QCMMopts.pkl")

	_path   = "/home/igorchem/VisMol/easyhybrid/pDynamoMethods/TIMtest_US1Dsimpledistance"
	potmean = PMF( proj.cSystem, _path, "TIM_FreeEnergy1DsimpleDistance" )
	potmean.CalculateWHAM(16,0,300.15)
	proj.FinishRun()
	
	
	_path   = "/home/igorchem/VisMol/easyhybrid/pDynamoMethods/TIMtest_US1DmultipleDistance"
	potmeanb = PMF( proj.cSystem, _path, "TIM_FreeEnergy1DmultipleDistance" )
	potmeanb.CalculateWHAM(16,0,300.15)
	potmeanb.FinishRun()
	
#=====================================================
def PotentialOfMeanField2D(self):
	'''
	'''
	proj=SimulationProject("TIMtest_PMF2D")		
	proj.LoadSystemFromSavedProject("TIMTest_QCMMopts.pkl")

	_path   = "/home/igorchem/VisMol/easyhybrid/pDynamoMethods/TIMtest_US2Dsimpledistance"
	potmean = PMF( proj.cSystem, _path, "TIM_FreeEnergy2DsimpleDistance" )
	potmean.CalculateWHAM(10,10,300.15)
	proj.FinishRun()
	
	
	_path   = "/home/igorchem/VisMol/easyhybrid/pDynamoMethods/TIMtest_US2DmixedDistance"
	potmeanb = PMF( proj.cSystem, _path, "TIM_FreeEnergy2DmixedDistance" )
	potmeanb.CalculateWHAM(10,10,300.15)
	potmeanb.FinishRun()

	_path   = "/home/igorchem/VisMol/easyhybrid/pDynamoMethods/TIMtest_US2DmultipleDistance"
	potmeanb = PMF( proj.cSystem, _path, "TIM_FreeEnergy2DmultipleDistance" )
	potmeanb.CalculateWHAM(10,10,300.15)
	potmeanb.FinishRun()

#=====================================================
def ReacCoordSearchers(self):
	pass 

#=====================================================
def SMOEnergyRef(self):
	'''
	'''
	'''
	proj=SimulationProject("TIMtest_SMOref")		
	proj.LoadSystemFromSavedProject("TIMTest_QCMMopts.pkl")
	methods= ""
	_path = os.path.join( os.getcwd(), "TIMtest_QCMM_Scans2DmultipleDistance","ScanTraj.ptGeo")
	parameters = { 'xnbins':6			,
				   "ynbins":6			,
				   "Scr_folder":_path   ,
				   "charge":-3		    ,
				   "multiplicity":1 	,
				   "methods_lists":methods,					   
				   "NmaxThreads":8 		,
				   "Software":"ORCA"	,
				   "rc1_min":-1.0       ,			
				   "rc2_min":-1.0       ,	
				   "rc2_max":-0.4		,
				   "rc1_max":-0.2       ,
				   "rc1_name":"Reaction Coord1",
				   "rc2_name":"Reaction Coord2"
				 }

	proj.RunSimulation(parameters,"Energy_Refinement")
	'''
	proj=SimulationProject("TIMtest_ORCAref")		
	proj.LoadSystemFromSavedProject("TIMTest_QCMMopts.pkl")
	methods= ""
	_path = os.path.join( os.getcwd(), "TIMtest_QCMM_Scans2DmultipleDistance","ScanTraj.ptGeo")
	parameters = { 'xnbins':6			,
				   "ynbins":6			,
				   "Scr_folder":_path   ,
				   "charge":-3		    ,
				   "multiplicity":1 	,
				   "orca_method":methods,
				   "basis":"6-31G*"		,					   
				   "NmaxThreads":8 		,
				   "Software":"ORCA"	,
				   "rc1_min":-1.0       ,			
				   "rc2_min":-1.0       ,	
				   "rc2_max":-0.4		,
				   "rc1_max":-0.2       ,
				   "rc1_name":"Reaction Coord1",
				   "rc2_name":"Reaction Coord2"
				 }

	proj.RunSimulation(parameters,"Energy_Refinement")

#=====================================================
def Thermodynamics(self):
	pass

#=====================================================
def GetOrbitalsInfo(self):
	pass


#============================================================================
if __name__ == "__main__":
	logFile.Header()

	logFile.Footer()
