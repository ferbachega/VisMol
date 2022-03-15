#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = Tests.py

#-------------------------------------------------------------
#Script file for the EasyHybrid 3.0 Core functionalities tests 

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#=================================================================
import os, glob, sys 

os.environ['MPLCONFIGDIR'] = '/tmp'


#path fo the core python files on your machine
sys.path.append("/home/igorchem/VisMol/easyhybrid/pDynamoMethods") 
#------------------------------------------------------
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

from EnergyAnalysis import EnergyAnalysis
from TrajectoryAnalysis import TrajectoryAnalysis
from ReactionCoordinate import *
#-------------------------------------------------------------------
#path for the required files on the examples folder of EasyHynrid 3.0
easyhybrid   = "/home/igorchem/VisMol/easyhybrid"
ex_path      = "/home/igorchem/VisMol/examples/"
scratch_path = os.path.join(easyhybrid,"TestsScratch")
timTop       = os.path.join(ex_path,"TIM","7tim.top")
timCrd       = os.path.join(ex_path,"TIM","7tim.crd")
balapkl      = os.path.join(ex_path,"bala","bAla.pkl")
#--------------------------------------------------------
if not os.path.exists(scratch_path):
	os.makedirs(scratch_path)

#*************************************************************************
def SetMMsytem():
	'''
	TESTED !!!
	Use the methods of the CoreInterface Class to set the Triosephosphate isomerase System
	'''
	#--------------------------------------------------------------------
	proj = SimulationProject( os.path.join(scratch_path,"MM_SetUp") )
	proj.LoadSystemFromForceField(timTop,timCrd)	
	#optimize full system
	parameters_a = {"maxIterations":1000,"rmsGradient":1}
	_plotParameters = None
	proj.RunSimulation(parameters_a,"Geometry_Optimization",_plotParameters)
	_pattern = "*:LIG.248:C02"
	#prune system to spherical selection
	proj.SphericalPruning(_pattern,25.0)
	#Fixed external atoms
	proj.SettingFixedAtoms(_pattern,20.0)
	parameters_b = {"maxIterations":1000,"rmsGradient":0.1}
	#otimize pruned systems
	proj.RunSimulation(parameters_b,"Geometry_Optimization",_plotParameters)
	#seve a pkl with the MM model defined for the pruned system 
	proj.SaveProject()
	proj.FinishRun() 
#=====================================================
def MMMD_Algorithms():
	'''
	TESTED !!!
	Molecular Mechanics Molecular Dynamics 
	'''
	#---------------------------------------------------
	#If the pkl with the Pruned MM system does not exist, generate it
	if not os.path.exists( os.path.join(scratch_path,"MM_SetUp.pkl") ):
		SetMMsytem()
	#------------------------------------------------
	proj=SimulationProject( os.path.join(scratch_path,"MM_MDAlgs") )	
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path, "MM_SetUp.pkl") )
	refcrd3 = Clone(proj.cSystem.coordinates3)
	#-----------------------------------------------
	integrators = ["Verlet", "LeapFrog", "Langevin"]	
	_plotParameters = {"show":True}
	#------------------------------------------------
	#loop to execute the available intgrators on pDynamo
	for integrator in integrators:
		#Non-optional parameters for molecular dynamics simulation preset
		parameters = {"protocol":"production"     ,
					  "production_nsteps":2000    ,
					  "equilibration_nsteps":1000 ,
					  "MD_method":integrator      }
		proj.RunSimulation(parameters,"Molecular_Dynamics",_plotParameters)
		proj.cSystem.coordinates3 = refcrd3
	#_--------------
	proj.FinishRun() 			
#=====================================================
def MMMD_Protocols():
	'''
	TESTED !!!
	Molecular Mechanics Molecular Dynamics 
	'''
	#-----------------------------------------------
	#If the pkl with the Pruned MM system does not exist, generate it
	if not os.path.exists( os.path.join( scratch_path,"MM_SetUp.pkl") ):
		SetMMsytem()
	#-----------------------------------------------
	proj=SimulationProject( os.path.join( scratch_path, "MM_MD_protocols") )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"MM_SetUp.pkl") )
	#-----------------------------------------------
	protocols = [ "heating", "equilibration", "production"]	
	_plotParameters = {"show":True}
	refcrd3 = Clone(proj.cSystem.coordinates3)
	#-------------------------
	for protocol in protocols:
		#Non-optional parameters for molecular dynamics simulation preset
		parameters = {"protocol":protocol         ,
					  "production_nsteps":2000    ,
					  "equilibration_nsteps":1000 ,
					  "MD_method":"LeapFrog"	  }
		proj.RunSimulation(parameters,"Molecular_Dynamics",_plotParameters)
		proj.cSystem.coordinates3 = refcrd3
	#---------------
	proj.FinishRun()			
#=====================================================
def QCMM_Energies():
	'''
	TESTED !!!
	Setting Quantum chemical model and saving the region.
	Test single point energy calculations using all Semiempirical methods provided in pDynamo.
	'''
	#-----------------------------------------------
	#If the pkl with the Pruned MM system does not exist, generate it
	if not os.path.exists( os.path.join( scratch_path, "MM_SetUp.pkl") ):
		SetMMsytem()
	#-----------------------------------------------	
	proj= SimulationProject( os.path.join( scratch_path,"QCMM_SMOs") )
	proj.LoadSystemFromSavedProject( os.path.join( scratch_path,"MM_SetUp.pkl") )
	#Defining QC region
	lig = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.248:*")
	glu = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:*")
	his = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:*")
	selections= [ lig, glu, his ]
	#-----------------------------
	#List of internal Hamiltonians provided by pDynamo with parameters for all atoms
	SMOmodels = ["am1","am1dphot","pddgpm3","pm3","rm1","pm6"]
	#-----------------------------
	#saving qc/mm setup
	for smo in SMOmodels:
		proj.SetSMOHybridModel( smo, selections, -3, 1 )
		proj.RunSinglePoint()	
	#saving qc/mm setup
	proj.SaveProject()
	proj.FinishRun()
#=====================================================
def QCMM_DFTBplus():
	'''
	TESTED !!!
	Setting quantum chemical system to Run interfaced with DFTB+
	'''
	#-----------------------------------------------
	#If the pkl with the Pruned QC system does not exist, generate it
	if not os.path.exists( os.path.join( scratch_path, "QCMM_SMOs.pkl") ):
		QCMM_Energies()
	#-----------------------------------------------
	proj= SimulationProject( os.path.join( scratch_path, "QCMM_DFTB") )
	proj.LoadSystemFromSavedProject( os.path.join( scratch_path,"QCMM_SMOs.pkl" ) )
	pureQCAtoms = list(proj.cSystem.qcState.pureQCAtoms)		
	proj.SetDFTBsystem(pureQCAtoms, -3, 1 )
	proj.FinishRun()
#=====================================================
def QCMM_Orca():
	'''
	TESTED !!!
	Test interface of ORCA
	'''
	#-----------------------------------------------
	#If the pkl with the Pruned QCsystem does not exist, generate it
	if not os.path.exists( os.path.join( scratch_path, "QCMM_SMOs.pkl") ):
		QCMM_Energies()
	#-----------------------------------------------
	proj= SimulationProject( os.path.join( scratch_path, "QCMM_ORCA") )
	proj.LoadSystemFromSavedProject( os.path.join( scratch_path,"QCMM_SMOs.pkl") )
	pureQCAtoms = list(proj.cSystem.qcState.pureQCAtoms)
	proj.SetOrcaSystem( "HF","6-31G*",pureQCAtoms,-3, 1 )
	proj.RunSinglePoint()
	proj.FinishRun()
#=====================================================
def QCMM_optimizations():
	'''
	TESTED !!!
	Test the Geometry optimizations algorithms with hybrid potential system
	'''
	#-----------------------------------------------
	#If the pkl with the Pruned QCsystem does not exist, generate it
	if not os.path.exists( os.path.join( scratch_path,"QCMM_SMOs.pkl") ):
		QCMM_Energies()
	#-----------------------------------------------
	proj= SimulationProject(  os.path.join( scratch_path,"QCMMopts") )
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMM_SMOs.pkl") )
	initialCoords = Clone(proj.cSystem.coordinates3)
	#Quasi-Newton muito demorado
	#---------------------------
	algs = ["ConjugatedGradient",
			"LFBGS"             ,
			"SteepestDescent"   ,
			"FIRE"              ]	
	#problems saving trajectory! But geometry opt working
	_plotParameters = None
	#---------------------------------------------------
	for alg in algs:
		parameters = {	"maxIterations":1000    ,
						"rmsGradient":0.1       ,
						"optmizer":alg			}						
		proj.RunSimulation(parameters,"Geometry_Optimization",_plotParameters)
		proj.cSystem.coordinates3 = initialCoords;
	#Save QCMM optimezed System	
	proj.SaveProject()
	proj.FinishRun()	
#=====================================================
def QCMM_MD():
	'''
	TESTED !!!
	Test QCMM hybrid molecular dynamics
	'''
	#-----------------------------------------------
	#If the pkl with the Pruned QCsystem does not exist, generate it
	if not os.path.exists( os.path.join( scratch_path, "QCMMopts.pkl") ):
		QCMM_optimizations()
	#------------------------------------------------
	_plotParameters = {"show":True}
	proj=SimulationProject( os.path.join(scratch_path,"QCMM_MDs") )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	#testing qcmm MD 
	#--------------------------------------------------------
	parameters = {	"protocol":"production"      ,
					"production_nsteps":10000    ,
					"equilibration_nsteps":5000  ,
					"MD_method":"LeapFrog"       ,
					"sampling_factor":400        }
	#--------------------------------------------------------			
	proj.RunSimulation(parameters,"Molecular_Dynamics",_plotParameters)
	proj.FinishRun()
#=====================================================
def QCMM_MDrestricted():
	'''
	TESTED!!!
	Test QCMM molecular dnamics with restricted reaction coorinates
	'''
	#---------------------------------------------
	#If the pkl with the Pruned QCsystem does not exist, generate it
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMM_optimizations()
	#---------------------------------------------
	proj=SimulationProject( os.path.join(scratch_path,"QCMM_restricted") )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )		
	#testing qcmm MD 
	#---------------------------------------------------------------
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")
	atomsf = [ atom1[0], atom2[0], atom3[0] ]
	atom6 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:O06")
	atom5 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")	
	atomss = [ atom4[0], atom5[0], atom6[0] ]
	#-----------------------------------------------------------------
	parameters = {	"protocol":"production"     ,
					"production_nsteps":1500    ,
					"equilibration_nsteps":500  ,
					"MD_method":"LeapFrog"      ,
				 	"atoms_M1":atomsf           ,
				 	"atoms_M2":atomss           ,
				 	"force_constant":100.0      ,
				 	"ndim":2                    ,
				 	"MC_RC1":True               ,
				 	"MC_RC2":True               ,
				 	"sampling_factor":100       }
	_plotParameters = {"show":True} 
	#------------------------------------------------
	#Run simulation	
	proj.RunSimulation(parameters,"Restricted_Molecular_Dynamics",_plotParameters)
	proj.FinishRun()		
#=====================================================
def QCMMScanSimpleDistance(_nsteps,_dincre,name="Default"):
	'''
	TESTED !
	Test QCMM one-dimensional reaction scans, using simple
	'''
	#---------------------------------------------
	#If the pkl with the Pruned QCsystem does not exist, generate it
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMM_optimizations()
	#---------------------------------------------
	_scanFolder = "QCMM_SCAN1D_simple_distance"
	if not name == "Default":
		_scanFolder = name
	proj=SimulationProject( os.path.join(scratch_path,_scanFolder) )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	_plotParameters = { "contour_lines":15 }
	#setting atoms for scan
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atomsf = [ atom1[0], atom2[0] ] 
	#--------------------------------------------
	#set parameters for relaxed surface scan
	parameters = { "ATOMS_RC1"  :atomsf		,
				   "dincre_RC1" :_dincre	,
				   "nSteps_RC1":_nsteps     ,
				   "ndim"      :1    		,
				   "MC_RC1"    :"true"      ,
				   "force_constant":4000.0	}
	
	#QCMM scans simple distances
	proj.RunSimulation(parameters,"Relaxed_Surface_Scan",_plotParameters)		
	proj.FinishRun()
#===============================================================
def QCMMScanMultipleDistance(_nsteps,_dincre,name="Default"):
	'''
	Test QCMM one-dimensional reaction scans, using multiple distance 
	'''
	#---------------------------------------------
	#If the pkl with the Pruned QCsystem does not exist, generate it
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMMoptimizations()
	#---------------------------------------------	
	_scanFolder = "QCMM_SCAN1D_multiple_distance"
	if not name == "Default":
		_scanFolder = name
	proj=SimulationProject( os.path.join(scratch_path,_scanFolder) )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	#---------------------------------------------
	#setting atoms for scan
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")
	atomsf = [ atom1[0], atom2[0], atom3[0] ] 
	#setting parameters
	_plotParameters = { "contour_lines":15 }
	parameters = { "ATOMS_RC1":atomsf     ,
				   "dincre_RC1":_dincre   ,
				   "nSteps_RC1":_nsteps   , 
				   "ndim":1               ,
				   "MC_RC1":"true"        ,
				   "force_constant":4000.0}
    #run the simulation
    #---------------------------------------------------------------------
	proj.RunSimulation(parameters,"Relaxed_Surface_Scan",_plotParameters)		
	proj.FinishRun()
#=====================================================
def Scan1D_Dihedral(_nsteps,name="Default"):
	'''
	Test relaxed scan with dihedral reaction coordinates
	'''
	#---------------------------------------------	
	_scanFolder = "SCAN1D_dihedral"
	if not name == "Default":
		_scanFolder = name
	proj=SimulationProject( os.path.join(scratch_path,_scanFolder) )		
	proj.LoadSystemFromSavedProject( balapkl )
	proj.cSystem.Summary()
	#---------------------------------------------
	#setting atoms for scan
	
	atomsf = [ 4, 6,  8, 14] 
	#setting parameters
	_plotParameters = { "contour_lines":15 }
	parameters = { "ATOMS_RC1":atomsf     ,
				   "nSteps_RC1":_nsteps   ,
				   "rc_type_1" :"dihedral", 
				   "ndim":1               ,
				   "MC_RC1":"true"        ,
				   "force_constant":25.0}
    #run the simulation
    #---------------------------------------------------------------------
	proj.RunSimulation(parameters,"Relaxed_Surface_Scan",_plotParameters)		
	proj.FinishRun()
#=====================================================
def Scan2D_Dihedral(_xnsteps,_ynsteps,name="Default"):
	'''
	Test relaxed scan with dihedral reaction coordinates
	'''
	#---------------------------------------------	
	_scanFolder = "SCAN2D_dihedral"
	if not name == "Default":
		_scanFolder = name
	proj=SimulationProject( os.path.join(scratch_path,_scanFolder) )		
	proj.LoadSystemFromSavedProject( balapkl )
	proj.cSystem.Summary()
	#---------------------------------------------
	#setting atoms for scan
	
	atomsf = [ 4, 6,  8, 14] 
	atomss = [ 6, 8, 14, 16]
	#setting parameters
	_plotParameters = { "contour_lines":10 }
	parameters = { "ATOMS_RC1":atomsf     ,
				   "ATOMS_RC2":atomss     ,
				   "nSteps_RC1":_xnsteps  ,
				   "nSteps_RC2":_ynsteps  ,
				   "rc_type_1" :"dihedral", 
				   "rc_type_2" :"dihedral", 
				   "ndim":2               ,
				   "force_constant_1":25.0,
				   "force_constant_2":25.0,
				   "NmaxThreads":4        }
    #run the simulation
    #---------------------------------------------------------------------
	proj.RunSimulation(parameters,"Relaxed_Surface_Scan",_plotParameters)		
	proj.FinishRun()
#=====================================================
def QCMMScan2DsimpleDistance(_xnsteps,_ynsteps,_dincrex,_dincrey,name="Default"):
	'''
	Evaluate tests for two-dimensional surface scans, with reaction coordinates set as simple distances
	'''	
	#---------------------------------------------
	#If the pkl with the Pruned QCsystem does not exist, generate it
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMMoptimizations()
	_scanFolder = "QCMM_Scan2D_simple_distance"
	if not name == "Default":
		_scanFolder = name
	proj=SimulationProject( os.path.join(scratch_path,_scanFolder) )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	#----------------------------------------
	#Set atoms for th reaction coordinates
	atom1  = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2  = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom5  = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4  = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")
	atomsf = [ atom1[0],atom2[0] ] 	
	atomss = [ atom4[0], atom5[0] ]
	#----------------------------------------
	_plotParameters = { "contour_lines":15 }
	parameters = { 'ATOMS_RC1':atomsf	  ,
				   'ATOMS_RC2':atomss	  ,
				   'dincre_RC1':_dincrex  ,
				   'dincre_RC2':_dincrey  , 
				   "nSteps_RC1":_xnsteps  ,
				   "nSteps_RC2":_ynsteps  , 
				   "ndim": 2 			  ,				 
				   "NmaxThreads":        8}

	proj.RunSimulation(parameters,"Relaxed_Surface_Scan",_plotParameters)		
	proj.FinishRun()	
#=====================================================
def QCMMScan2DmixedDistance(_xnsteps,_ynsteps,_dincrex,_dincrey,name="Default"):
	'''
	'''
	#-------------------------------------------------
	#If the pkl with the Pruned QCsystem does not exist, generate it
	if not os.path.exists( os.path.join(scratch_path, "QCMMopts.pkl") ):
		QCMMoptimizations()
	_scanFolder = "QCMM_Scan2D_mixed_distance"
	if not name == "Default":
		_scanFolder = name
	proj=SimulationProject( os.path.join(scratch_path, _scanFolder ) )	
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	#--------------------------------------------------
	#set atoms for reaction coordinates
	atom1  = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2  = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3  = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")	
	atom5  = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4  = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")
	atomsf = [ atom1[0],atom2[0], atom3[0] ] 	
	atomss = [ atom4[0], atom5[0] ]
	#--------------------------------------------------
	_plotParameters = { "contour_lines":15 }

	parameters = { "ATOMS_RC1":atomsf	  ,
				   "ATOMS_RC2":atomss	  ,
				   "dincre_RC1":_dincrex  ,
				   "dincre_RC2":_dincrey  , 
				   "nSteps_RC1":_xnsteps  ,
				   "nSteps_RC2":_ynsteps  , 
				   "ndim": 2 			  ,
				   "MC_RC1":		"true",
				   "NmaxThreads":        8,
				 }
	#--------------------------------------------------
	proj.RunSimulation(parameters,"Relaxed_Surface_Scan",_plotParameters)		
	proj.FinishRun()
#=====================================================
def QCMMScan2DmultipleDistance(_xnsteps,_ynsteps,_dincrex,_dincrey,name="Default"):
	'''
	'''
	#-------------------------------------------------
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMMoptimizations()
	_scanFolder = "QCMM_Scan2D_multiple_distance"
	if not name == "Default":
		_scanFolder = name
	proj=SimulationProject( os.path.join(scratch_path,_scanFolder) )
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	#--------------------------------------------------

	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")	
	atom6 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:O06")
	atom5 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")
	#--------------------------------------------------
	atomsf = [ atom1[0], atom2[0], atom3[0] ] 
	atomss = [ atom4[0], atom5[0], atom6[0] ]
	#--------------------------------------------------
	_plotParameters = { "contour_lines":15 }
	#--------------------------------------------------
	parameters = { 'ATOMS_RC1':atomsf	  ,
				   'ATOMS_RC2':atomss	  ,
				   'dincre_RC1':_dincrex  ,
				   'dincre_RC2':_dincrey  , 
				   "nSteps_RC1":_xnsteps  ,
				   "nSteps_RC2":_ynsteps  , 
				   "ndim": 2 			  ,
				   "MC_RC1":		"true",
				   "MC_RC2":		"true",
				   "NmaxThreads":        8}
				 
	proj.RunSimulation(parameters,"Relaxed_Surface_Scan",_plotParameters)		
	proj.FinishRun()	
#=====================================================
def QCMMScans2D_Adaptative(_xnsteps,_ynsteps,_dincrex,_dincrey):
	'''
	'''
	#-------------------------------------------------
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMMoptimizations()
	proj=SimulationProject( os.path.join(scratch_path,"QCMM_Scan2D_adaptative")	)	
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	#--------------------------------------------------
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")	
	atom6 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:O06")
	atom5 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")
	#--------------------------------------------------
	atomsf = [ atom1[0], atom2[0], atom3[0] ] 
	atomss = [ atom4[0], atom5[0], atom6[0] ]
	#--------------------------------------------------
	_plotParameters = { "contour_lines":15,"show":"true"}
	#--------------------------------------------------
	parameters = { 'ATOMS_RC1':atomsf	  ,
				   'ATOMS_RC2':atomss	  ,
				   'dincre_RC1':_dincrex  ,
				   'dincre_RC2':_dincrey  , 
				   "nSteps_RC1":_xnsteps  ,
				   "nSteps_RC2":_ynsteps  , 
				   "ndim": 2 			  ,
				   "MC_RC1":		"true",
				   "MC_RC2":		"true",
				   "NmaxThreads":        8,
				   "adaptative":"true"    }
	proj.RunSimulation(parameters,"Relaxed_Surface_Scan",_plotParameters)		
	proj.FinishRun()	
#=====================================================
def FreeEnergy1DSimpleDistance(nsteps):
	'''
	Test simulations to get the Free energy of one-dimensional reaction scan using umbrella sampling and WHAM
	'''
	#-------------------------------------------------
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMM_optimizations()
	proj=SimulationProject( os.path.join(scratch_path,"FE_simple_distance") )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	#-------------------------------------------------
	_name = "SCAN1D_4FEcalculations_simple_distance"
	_path = os.path.join( os.path.join(scratch_path,_name,"ScanTraj.ptGeo") )
	QCMMScanSimpleDistance(10,0.2,name=_name)
	#-------------------------------------------------
	atom1  = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2  = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atomsf = [ atom1[0], atom2[0] ]		
	_plotParameters = { "contour_lines":15,"xwindows":10,"ywindows":0,"crd1_label":"Reaction Coordinate"}
	#-------------------------------------------------	
	parameters = { "ATOMS_RC1":atomsf			  ,
				   "ndim": 1 					  ,
				   "sampling_factor":nsteps/10	  ,
				   "equilibration_nsteps":nsteps/2,
				   "production_nsteps":nsteps	  ,
				   "source_folder":_path 		  ,
				   "MD_method":"LeapFrog"		  ,
				   "MC_RC1":"true"				  ,
				   "NmaxThreads":8 				  }
	#-------------------------------------------------
	#RUN umbrella sampling
	proj.RunSimulation(parameters,"Umbrella_Sampling",_plotParameters)	
	#-------------------------------------------------
	#path for the ptRes files
	_path = os.path.join( scratch_path, "FE_simple_distance" )	
	parameters = { "source_folder":_path ,
				   "xnbins":20           ,
				   "ynbins":0            ,
				   "temperature":300.15	 }
	#RUN WHAM, calculate PMF and free energy
	proj.RunSimulation(parameters,"PMF_Analysis",_plotParameters)
#=================================================================
def FreeEnergy1DMultipleDistance(nsteps):
	'''
	Test simulations to get the Free energy of one-dimensional reaction scan using umbrella sampling and WHAM
	'''
	#-----------------------------------------------------------
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMM_optimizations()
	#-----------------------------------------------------------
	proj=SimulationProject( os.path.join(scratch_path,"FE_multiple_distance") )	
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	#-------------------------------------------------
	_name = "SCAN1D_4FEcalculations_multiple_distance"
	_path = os.path.join( os.path.join(scratch_path,_name,"ScanTraj.ptGeo" ) )
	QCMMScanMultipleDistance(10,0.2,name=_name)	
	#-------------------------------------------------
	#set atoms for the reaction coordinates
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")
	#-------------------------------------------------
	_plotParameters = { "contour_lines":15,"xwindows":10,"ywindows":0,"crd1_label":"Reaction Coordinate"}
	atomsf = [ atom1[0], atom2[0], atom3[0] ] 
	#-------------------------------------------------	
	parameters = { 'ATOMS_RC1':atomsf				,
				   "ndim": 1 						,
				   "sampling_factor":nsteps/10    	,
				   "equilibration_nsteps":nsteps/2 	,
				   "production_nsteps":nsteps		,
				   "source_folder":_path 	    	,
				   "MD_method":"LeapFrog"			,
				   "MC_RC1":"true"					,
				   "NmaxThreads":8 				    }
	#-------------------------------------------------
	#Run umbrella sampling 
	proj.RunSimulation(parameters,"Umbrella_Sampling",_plotParameters)

	_path = os.path.join( scratch_path,"FE_multiple_distance")
	parameters = { "source_folder":_path,
				   "xnbins":20          ,
				   "ynbins":0           ,
				   "temperature":300.15 }
	#-------------------------------------------------
	#Run umbrella sampling 
	proj.RunSimulation(parameters,"PMF_Analysis",_plotParameters)
	proj.FinishRun()
#=====================================================
def UmbrellaSampling1Drestart(nsteps):
	'''
	Test the restart utility of the Umbrella sapling preset simulation
	'''
	#-----------------------------------------------------------
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMM_optimizations()
	#-----------------------------------------------------------
	proj=SimulationProject( os.path.join(scratch_path,"UmbrellaSampling_Restart") )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )

	_name = "SCAN1D_4US_restart_test"
	QCMMScanMultipleDistance(6,0.2,name=_name)
	_path = os.path.join( scratch_path,_name,"ScanTraj.ptGeo" )
	
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")
	#-------------------------------------------------
	_plotParameters = { "contour_lines":15}
	atomsf = [ atom1[0], atom2[0], atom3[0] ] 
	#-------------------------------------------------	
	parameters = { 'ATOMS_RC1':atomsf				,
				   "ndim": 1 						,
				   "sampling_factor":nsteps/10    	,
				   "equilibration_nsteps":nsteps/2 	,
				   "production_nsteps":nsteps		,
				   "source_folder":_path 	    	,
				   "MD_method":"LeapFrog"			,
				   "MC_RC1":"true"					,
				   "restart":"true"					,
				   "NmaxThreads":8 				    }
	#-------------------------------------------------
	#Run first umbrella sampling 
	proj.RunSimulation(parameters,"Umbrella_Sampling",_plotParameters)

	_pathpmf = os.path.join(scratch_path,"UmbrellaSampling_Restart")
	_PMFparameters = { "source_folder":_pathpmf,
				   "xnbins":10          ,
				   "ynbins":0           ,
				   "temperature":300.15 }
	#-------------------------------------------------
	#Run umbrella sampling 
	proj.RunSimulation(_PMFparameters,"PMF_Analysis",_plotParameters)
	#*************************************************************
	QCMMScanMultipleDistance(10,0.2)
	parameters = { 'ATOMS_RC1':atomsf				,
				   "ndim": 1 						,
				   "sampling_factor":nsteps/10    	,
				   "equilibration_nsteps":nsteps/2 	,
				   "production_nsteps":nsteps		,
				   "source_folder":_path 	    	,
				   "MD_method":"LeapFrog"			,
				   "MC_RC1":"true"					,
				   "restart":"true"                 ,
				   "NmaxThreads":8 				    }

	proj.RunSimulation(parameters,"Umbrella_Sampling",_plotParameters)

	proj.RunSimulation(_PMFparameters,"PMF_Analysis",_plotParameters)
	proj.FinishRun()
#=====================================================
def FreeEnergy2DsimpleDistance(nsteps):
	'''
	'''
	proj=SimulationProject( os.path.join(scratch_path,"FE_2D_simple_distance") )
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMM_optimizations()		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	#-------------------------------------------------------------------
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")	
	atom6 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:O06")
	atom5 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")
	atomsf = [ atom1[0], atom2[0] ] 
	atomss = [ atom4[0], atom5[0] ]
	#------------------------------------------------------------------
	_name = "SCAN2D_4FEcalculations_simple_distance"
	_path = os.path.join( scratch_path,_name,"ScanTraj.ptGeo")
	if not os.path.exists(_path):
		QCMMScan2DsimpleDistance(6,6,0.2,0.2,name=_name)	
	
	_plotParameters = { "contour_lines":12,"xwindows":6,"ywindows":6,"crd1_label":"Reaction Coordinate"}	

	parameters = { "ATOMS_RC1":atomsf				,
				   "ATOMS_RC2":atomss				,
				   "ndim": 2 						,
				   "sampling_factor":nsteps/10    	,
				   "equilibration_nsteps":nsteps/2  ,
				   "production_nsteps":nsteps		,
				   "source_folder":_path 			,
				   "MD_method":"LeapFrog"			,
				   "MC_RC1":"true"					,
				   "MC_RC2":"true"					,
				   "restart":"True"                 ,
				   "NmaxThreads":8 					}
	
	proj.RunSimulation(parameters,"Umbrella_Sampling",_plotParameters)
	_path = os.path.join( scratch_path, "FE_2D_simple_distance")	
	parameters = { "source_folder":_path ,
				   "xnbins":10           ,
				   "ynbins":10           ,
				   "temperature":300.15	 }
	#RUN WHAM, calculate PMF and free energy
	proj.RunSimulation(parameters,"PMF_Analysis",_plotParameters)
	proj,FinishRun()
#=====================================================
def FreeEnergy2DmixedDistance(nsteps):
	'''
	'''
	proj=SimulationProject( os.path.join(scratch_path,"FE_2D_mixed_distance") )	
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMM_optimizations()
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")	
	atom6 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:O06")
	atom5 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")
	atomsf = [ atom1[0], atom2[0], atom3[0] ] 
	atomss = [ atom4[0], atom5[0] ]

	_name = "SCAN2D_4FEcalculations_mixed_distance"
	_path = os.path.join( scratch_path,_name,"ScanTraj.ptGeo")
	if not os.path.exists(_path):
		QCMMScan2DmixedDistance(6,6,0.2,0.2,name=_name)

	_plotParameters = { "contour_lines":15,"xwindows":6,"ywindows":6,"crd1_label":"Reaction Coordinate"}

	parameters = { 'ATOMS_RC1':atomsf			,
				   "ATOMS_RC2":atomss			,
				   "ndim": 2 					,
				   "sampling_factor":nsteps/10 	,
				   "equilibration_nsteps":nsteps/2,
				   "production_nsteps":nsteps	,
				   "source_folder":_path 		,
				   "MD_method":"LeapFrog"		,
				   "MC_RC1":"true"				,
				   "MC_RC2":"true"				,
				   "restart":"true"             , 
				   "NmaxThreads":8 				}
	
	proj.RunSimulation(parameters,"Umbrella_Sampling",_plotParameters)

	_path = os.path.join( scratch_path, "FE_2D_mixed_distance")	
	parameters = { "source_folder":_path ,
				   "xnbins":10           ,
				   "ynbins":10           ,
				   "temperature":300.15	 }
	#RUN WHAM, calculate PMF and free energy
	proj.RunSimulation(parameters,"PMF_Analysis",_plotParameters)
	proj.FinishRun()
#=====================================================
def FreeEnergy2DmultipleDistance(nsteps):
	'''
	'''
	proj=SimulationProject( os.path.join(scratch_path,"FE_2D_multiple_distance") )	
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMM_optimizations()
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")
	
	atom6 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:O06")
	atom5 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")
	atomsf = [ atom1[0], atom2[0], atom3[0] ] 
	atomss = [ atom4[0], atom5[0], atom6[0] ]

	_plotParameters = { "contour_lines":15,"xwindows":6,"ywindows":6,"crd1_label":"Reaction Coordinate"}
	
	_name = "SCAN2D_4FEcalculations_multiple_distance"
	_path = os.path.join( scratch_path,_name,"ScanTraj.ptGeo")
	if not os.path.exists(_path):
		QCMMScan2DmultipleDistance(6,6,0.2,0.2,name=_name)

	parameters = { "ATOMS_RC1":atomsf				,
				   "ATOMS_RC2":atomss				,
				   "ndim": 2 						,
				   "sampling_factor":nsteps/10		,
				   "equilibration_nsteps":nsteps/2 	,
				   "production_nsteps":nsteps		,
				   "source_folder":_path 			,
				   "MD_method":"LeapFrog"			,
				   "MC_RC1":"true"					,
				   "MC_RC2":"true"					,
				   "NmaxThreads":8 	     			}
	
	proj.RunSimulation(parameters,"Umbrella_Sampling",_plotParameters)

	_path = os.path.join( scratch_path, "FE_2D_multiple_distance")	
	parameters = { "source_folder":_path ,
				   "xnbins":12           ,
				   "ynbins":12           ,
				   "temperature":300.15	 }
	#RUN WHAM, calculate PMF and free energy
	proj.RunSimulation(parameters,"PMF_Analysis",_plotParameters)
	proj.FinishRun()
#=====================================================
def EnergyAnalysisPlots():
	'''
	General tests of plots from energies analysis from our simulations
	'''
	#---------------------------------------------------------------------------------
	#LOAD system
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMM_optimizations()
	proj=SimulationProject( os.path.join(scratch_path,"Energy Plots") )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	#----------------------------------------------------------------
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")	
	atom6 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:O06")
	atom5 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")
	#setting reaction coordinates for ploting labels
	a1 = [atom1[0],atom2[0]]
	rc1_sd = ReactionCoordinate(a1,False)
	rc1_sd.SetInformation(proj.cSystem,0)
	a1 = [atom1[0],atom2[0],atom3[0]]
	rc1_md = ReactionCoordinate(a1,False)
	rc1_md.SetInformation(proj.cSystem,0)
	a2 = [atom4[0],atom5[0]]
	rc2_sd = ReactionCoordinate(a2,False)
	rc2_sd.SetInformation(proj.cSystem,0)
	a2 = [atom4[0],atom5[0],atom6[0]]
	rc2_md = ReactionCoordinate(a2,False)
	rc2_md.SetInformation(proj.cSystem,0)

	#----------------------------------------------------------------------------------
	#1D energy plot test 
	if not os.path.exists( os.path.join(scratch_path,"QCMM_SCAN1D_simple_distance") ):
		QCMMScanSimpleDistance(30,0.05)
	log_path = os.path.join(scratch_path,"QCMM_SCAN1D_simple_distance_energy.log")
	parameters      = {"xsize":30,"type":"1D","log_name":log_path}
	_plotParameters = {"show":True,"crd1_label":rc1_sd.label}
	proj.RunSimulation(parameters,"Energy_Plots",_plotParameters)
	#===================================================================================
	#2D PES plots 
	if not os.path.exists( os.path.join(scratch_path,"QCMM_Scan2D_simple_distance") ):
		QCMMScan2DsimpleDistance(12,12,0.15,0.15)
	log_path = os.path.join(scratch_path,"QCMM_Scan2D_simple_distance_energy.log")
	parameters      = {"xsize":12,"ysize":12,"type":"2D","log_name":log_path}
	_plotParameters = {"show":True,"crd1_label":rc1_sd.label,"crd2_label":rc2_sd.label,"contour_lines":12}
	proj.RunSimulation(parameters,"Energy_Plots",_plotParameters)
	#------------------------------------------------------------
	if not os.path.exists( os.path.join(scratch_path,"QCMM_Scan2D_mixed_distance") ):
		QCMMScan2DmixedDistance(12,12,0.15,0.15)
	log_path = os.path.join(scratch_path,"QCMM_Scan2D_mixed_distance_energy.log")
	parameters      = {"xsize":12,"ysize":12,"type":"2D","log_name":log_path}
	_plotParameters = {"show":True,"crd1_label":rc1_md.label,"crd2_label":rc2_sd.label,"contour_lines":12}
	proj.RunSimulation(parameters,"Energy_Plots",_plotParameters)
	#------------------------------------------------------------
	if not os.path.exists( os.path.join(scratch_path,"QCMM_Scan2D_multiple_distance") ):
		QCMMScan2DmultipleDistance(12,12,0.15,0.15)
	log_path = os.path.join(scratch_path,"QCMM_Scan2D_multiple_distance_energy.log")
	parameters      = {"xsize":12,"ysize":12,"type":"2D","log_name":log_path}
	_plotParameters = {"show":True,"crd1_label":rc1_md.label,"crd2_label":rc2_md.label,"contour_lines":12}
	proj.RunSimulation(parameters,"Energy_Plots",_plotParameters)
	#===================================================================================
	# 1D Free energy and PMF plots 
	if not os.path.exists( os.path.join(scratch_path,"FE_simple_distance") ):
		FreeEnergy1DSimpleDistance(600)
	log_pathFE  = os.path.join(scratch_path,"FE_simple_distance_FE.log")
	log_pathPMF = os.path.join(scratch_path,"FE_simple_distance.dat")
	parameters      = {"xsize":10,"type":"FE1D","log_name":log_pathFE} # xsize is for windowns size
	_plotParameters = {"show":True,"crd1_label":rc1_sd.label}
	proj.RunSimulation(parameters,"Energy_Plots",_plotParameters)
	parameters      = {"xsize":20,"type":"WHAM1D","log_name":log_pathPMF} # xsize is for bins size for PMF
	_plotParameters = {"show":True,"crd1_label":rc1_sd.label}
	proj.RunSimulation(parameters,"Energy_Plots",_plotParameters)
	#------------------------------------------------------------------------------------
	if not os.path.exists( os.path.join(scratch_path,"FE_multiple_distance") ):
		FreeEnergy1DMultipleDistance(600)
	log_pathFE  = os.path.join(scratch_path,"FE_multiple_distance_FE.log")
	log_pathPMF = os.path.join(scratch_path,"FE_multiple_distance.dat")
	parameters      = {"xsize":10,"type":"FE1D","log_name":log_pathFE} # xsize is for windowns size
	_plotParameters = {"show":True,"crd1_label":rc1_md.label}
	proj.RunSimulation(parameters,"Energy_Plots",_plotParameters)
	parameters      = {"xsize":20,"type":"WHAM1D","log_name":log_pathPMF} # xsize is for bins size for PMF
	_plotParameters = {"show":True,"crd1_label":rc1_md.label}
	proj.RunSimulation(parameters,"Energy_Plots",_plotParameters)
	#------------------------------------------------------------------------------------
	#====================================================================================
	# 2D Free Energy FES and PMF 
	if not os.path.exists( os.path.join(scratch_path,"FE_2D_simple_distance") ):
		FreeEnergy2DsimpleDistance(1000)
	log_pathFE  = os.path.join(scratch_path,"FE_2D_simple_distance_FE.log")
	log_pathPMF = os.path.join(scratch_path,"FE_2D_simple_distance.dat")
	parameters      = {"xsize":6,"ysize":6,"type":"FE2D","log_name":log_pathFE} # xsize is for windowns size
	_plotParameters = {"show":True,"crd1_label":rc1_sd.label,"crd2_label":rc2_sd.label,"contour_lines":12}
	proj.RunSimulation(parameters,"Energy_Plots",_plotParameters)
	parameters      = {"xsize":12,"ysize":12,"type":"WHAM2D","log_name":log_pathPMF} # xsize is for bins size for PMF
	_plotParameters = {"show":True,"crd1_label":rc1_sd.label,"crd2_label":rc2_sd.label,"contour_lines":12}
	proj.RunSimulation(parameters,"Energy_Plots",_plotParameters)
	#-------------------------------------------------------------------------------
	if not os.path.exists( os.path.join(scratch_path,"FE_2D_mixed_distance") ):
		FreeEnergy2DmixedDistance(1000)
	log_pathFE  = os.path.join(scratch_path,"FE_2D_mixed_distance_FE.log")
	log_pathPMF = os.path.join(scratch_path,"FE_2D_mixed_distance.dat")
	parameters      = {"xsize":6,"ysize":6,"type":"FE2D","log_name":log_pathFE} # xsize is for windowns size
	_plotParameters = {"show":True,"crd1_label":rc1_sd.label,"crd2_label":rc2_sd.label,"contour_lines":12}
	proj.RunSimulation(parameters,"Energy_Plots",_plotParameters)
	parameters      = {"xsize":12,"ysize":12,"type":"WHAM2D","log_name":log_pathPMF} # xsize is for bins size for PMF
	_plotParameters = {"show":True,"crd1_label":rc1_sd.label,"crd2_label":rc2_sd.label,"contour_lines":12}
	proj.RunSimulation(parameters,"Energy_Plots",_plotParameters)
	#-------------------------------------------------------------------------------
	if not os.path.exists( os.path.join(scratch_path,"FE_2D_multiple_distance") ):
		FreeEnergy2DmultipleDistance(1000)
	log_pathFE = os.path.join(scratch_path,"FE_2D_multiple_distance_FE.log")
	log_pathPMF = os.path.join(scratch_path,"FE_2D_multiple_distance.dat")
	parameters      = {"xsize":6,"ysize":6,"type":"FE2D","log_name":log_pathFE} # xsize is for windowns size
	_plotParameters = {"show":True,"crd1_label":rc1_sd.label,"crd2_label":rc2_sd.label,"contour_lines":12}
	proj.RunSimulation(parameters,"Energy_Plots",_plotParameters)
	parameters      = {"xsize":12,"ysize":12,"type":"WHAM2D","log_name":log_pathPMF} # xsize is for bins size for PMF
	_plotParameters = {"show":True,"crd1_label":rc1_sd.label,"crd2_label":rc2_sd.label,"contour_lines":12}
	proj.RunSimulation(parameters,"Energy_Plots",_plotParameters)
	#===========================================================================
	#internal refinement
#===============================================================================
def TrajectoryAnalysisPlots():
	'''
	Test the analysis and plotting for molecular dynamics simulations.
	'''
	# Root mean square and radius of gyration plots 
	# Distance analysis in restricted molecular dynamics 
	# Extraction of most representative frame based on the metrics: rmsd, rg and reaction coordinate distances
#===============================================================================
def ReacCoordSearchers():	
	'''
	Test simulation protocols to determine reaction paths
	'''
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMM_optimizations()
	proj=SimulationProject( os.path.join(scratch_path,"ReactionPathsSearchers") )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )

	#generate initial and final coordinates for NEB 
	#generate trajectory for SAW
		
	_name = "SCAN1D_4NEB_and_SAW"
	_path = os.path.join( os.path.join(scratch_path,_name,"ScanTraj.ptGeo") )
	if not os.path.exists(_path):
		QCMMScanMultipleDistance(30,0.05,name=_name)
	init_path = os.path.join( _path, "frame0.pkl")
	final_path = os.path.join( _path, "frame29.pkl")

	parameters      = {"init_coord":init_path,"final_coord":final_path,"NEB_nbins":16,"RMS_growing_intial_string":0.3}
	_plotParameters = None 

	proj.RunSimulation(parameters,"NEB",_plotParameters)
#================================================================================
def pDynamoEnergyRef_1D():
	'''
	Tested
	'''	
	proj=SimulationProject( os.path.join(scratch_path, "pDynamoSMO") )
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMM_optimizations()		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )

	methods = ["am1","rm1","pm3","pm6","pddgpm3"]

	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")	
	#setting reaction coordinates for ploting labels
	a1 = [ atom1[0],atom2[0],atom3[0] ]
	rc1_md = ReactionCoordinate(a1,False)
	rc1_md.SetInformation(proj.cSystem,0)	
	a2 = [atom4[0],atom5[0],atom6[0]]
	rc2_md = ReactionCoordinate(a2,False)
	rc2_md.SetInformation(proj.cSystem,0)


	_name = "SCAN1D_4Refinement"
	_path = os.path.join( os.path.join(scratch_path,_name,"ScanTraj.ptGeo") )
	if not os.path.exists(_path):
		QCMMScanMultipleDistance(30,0.05,name=_name)

	_plotParameters = {	"show":True,
						"crd1_label":rc1_md.label,
						"crd2_label":rc2_md.label,
						"contour_lines":12 ,
						"xlim_list": [-1.2,2.0] }

	parameters = { "xnbins":30			,
				   "ynbins":0			,
				   "source_folder":_path,
				   "out_folder":os.path.join(scratch_path, "SMO_EnergyRefinement"),
				   "charge":-3		    ,
				   "multiplicity":1 	,
				   "methods_lists":methods,					   
				   "NmaxThreads":8 		,
				   "Software":"pDynamo"	}

	proj.RunSimulation(parameters,"Energy_Refinement",_plotParameters)	
#=====================================================
def pDynamoEnergyRef_2D():
	'''

	'''
	proj=SimulationProject( os.path.join(scratch_path, "pDynamoSMO") )
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMM_optimizations()		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )

	methods = ["am1","rm1"]

	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")	
	atom6 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:O06")
	atom5 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")
	#setting reaction coordinates for ploting labels
	a1 = [atom1[0],atom2[0],atom3[0]]
	rc1_md = ReactionCoordinate(a1,False)
	rc1_md.SetInformation(proj.cSystem,0)	
	a2 = [atom4[0],atom5[0],atom6[0]]
	rc2_md = ReactionCoordinate(a2,False)
	rc2_md.SetInformation(proj.cSystem,0)

	_name = "SCAN2D_4Refinement"
	_path = os.path.join( os.path.join(scratch_path,_name,"ScanTraj.ptGeo") )
	if not os.path.exists(_path):
		QCMMScan2DmultipleDistance(6,6,0.1,0.1,name=_name)

	_plotParameters = {	"show":True              ,
						"crd1_label":rc1_md.label,
						"crd2_label":rc2_md.label,
						"contour_lines":12       ,
						"xlim_list": [-1.2,-0.3] ,
						"ylim_list": [-0.9,-0.2] }

	parameters = { "xnbins":6			,
				   "ynbins":6			,
				   "source_folder":_path,
				   "out_folder":os.path.join(scratch_path, "SMO2D_EnergyRefinement"),
				   "charge":-3		    ,
				   "multiplicity":1 	,
				   "methods_lists":methods,					   
				   "NmaxThreads":4		,
				   "Software":"pDynamo"	}

	proj.RunSimulation(parameters,"Energy_Refinement",_plotParameters)
#=====================================================
def DFTBplusEnergy():
	'''
	'''
	pass
#=====================================================
def MopacEnergyRef():
	'''
	'''
	proj=SimulationProject( os.path.join(scratch_path, "mopacSMO") )
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMM_optimizations()		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	methods = ["am1","pm3","pm6","pm7","rm1"]
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")	
	atom6 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:O06")
	atom5 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4 = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")
	#setting reaction coordinates for ploting labels
	a1 = [atom1[0],atom2[0],atom3[0]]
	rc1_md = ReactionCoordinate(a1,False)
	rc1_md.SetInformation(proj.cSystem,0)	
	a2 = [atom4[0],atom5[0],atom6[0]]
	rc2_md = ReactionCoordinate(a2,False)
	rc2_md.SetInformation(proj.cSystem,0)
	_name = "SCAN1D_4MopacRefinement"
	_path = os.path.join( os.path.join(scratch_path,_name,"ScanTraj.ptGeo") )
	
	if not os.path.exists(_path):
		QCMMScanMultipleDistance(20,0.08,name=_name)

	_plotParameters = {	"show":True,
						"crd1_label":rc1_md.label,
						"crd2_label":rc2_md.label,
						"contour_lines":12 ,
						"xlim_list": [-1.2,2.0] }

	parameters = { "xnbins":20			,
				   "ynbins":0			,
				   "mopac_keywords":["grad qmmm","ITRY=5000"] ,
				   "source_folder":_path,
				   "out_folder":os.path.join(scratch_path, "MOPAC_EnergyRefinement"),
				   "charge":-3		    ,
				   "multiplicity":1 	,
				   "methods_lists":methods,					   
				   "NmaxThreads":1 		,
				   "Software":"mopac"	}

	proj.RunSimulation(parameters,"Energy_Refinement",_plotParameters)	
	
#=====================================================
def ORCAEnergy():
	pass
#=====================================================
def Thermodynamics():
	pass
#=====================================================
if __name__ == "__main__":
	#MMMD_Algorithms()                         			#TESTED
	#MMMD_Protocols()									#TESTED
	#QCMM_Energies()									#TESTED
	#QCMM_DFTBplus()									#TESTED
	#QCMM_Orca()										#TESTED
	#QCMM_optimizations()								#TESTED
	#QCMM_MD()											#TESTED
	#QCMM_MDrestricted()								#TESTED
	#QCMMScanSimpleDistance(30,0.05)					#TESTED
	#QCMMScanMultipleDistance(30,0.05)					#TESTED
	#QCMMScan2DsimpleDistance(12,12,0.1,0.1)			#TESTED
	#QCMMScan2DmixedDistance(12,12,0.1,0.1)				#TESTED
	#QCMMScan2DmultipleDistance(12,12,0.1,0.1)			#TESTED
	#QCMMScans2D_Adaptative(12,12,0.2,0.2)				#TESTED
	#FreeEnergy1DSimpleDistance(600)					#TESTED
	#FreeEnergy1DMultipleDistance(600)					#TESTED
	#UmbrellaSampling1Drestart(500)						#TESTED
	#FreeEnergy2DsimpleDistance(500)					#TESTED
	#FreeEnergy2DmixedDistance(500)						#TESTED
	#FreeEnergy2DmultipleDistance(500)					#TESTED
	#pDynamoEnergyRef_1D()								#TESTED
	EnergyAnalysisPlots()								#TESTED
	#ReacCoordSearchers()								#NEB TESTED
	#MopacEnergyRef()									#TESTED
	#pDynamoEnergyRef_2D()								#TESTED
	#Scan1D_Dihedral(36)								#TESTED
	#Scan2D_Dihedral(12,12)