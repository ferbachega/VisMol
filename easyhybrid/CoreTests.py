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
VISMOL_HOME = os.environ.get('VISMOL_HOME')
#path fo the core python files on your machine
#para funcionar nas minhas máquinas, depois ver forma melhor de fazer
if not VISMOL_HOME == None:
	sys.path.append(os.path.join(VISMOL_HOME,"easyhybrid/pDynamoMethods") ) 
else:
	VISMOL_HOME = "/home/igorchem/VisMol"
	sys.path.append(os.path.join("/home/igorchem/VisMol/easyhybrid/pDynamoMethods") ) 
#------------------------------------------------------
from pBabel                    import *                                     
from pCore                     import *                                     
from pMolecule                 import *                    
from CoreInterface 			   import SimulationProject
from ReactionCoordinate import *
#-------------------------------------------------------------------
#path for the required files on the examples folder of EasyHynrid 3.0
easyhybrid   = os.path.join(VISMOL_HOME, "easyhybrid")
ex_path      = os.path.join(VISMOL_HOME, "examples/")

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
	parameters_a = {"simulation_type":"Geometry_Optimization","maxIterations":1000,"rmsGradient":1 }	
	proj.RunSimulation(parameters_a)
	#prune system to spherical selection
	_pattern = "*:LIG.248:C02"
	proj.SphericalPruning(_pattern,25.0)
	proj.SettingFixedAtoms(_pattern,20.0)
	parameters_b = {"simulation_type":"Geometry_Optimization","maxIterations":1000,"rmsGradient":1,"save_pdb":True}
	#otimize pruned systems
	proj.RunSimulation(parameters_b)
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
	parameters = {"protocol":"sampling","temperature": 315.15,"simulation_type":"Molecular_Dynamics"}
	#------------------------------------------------
	#loop to execute the available intgrators on pDynamo
	for integrator in integrators:
		parameters["MD_method"]       = integrator
		parameters["nsteps"]          = 1000
		parameters["sampling_factor"] = 0 
		proj.RunSimulation(parameters)
		parameters["MD_method"]       = integrator
		parameters["nsteps"]   		  = 2000
		parameters["sampling_factor"] = 50 
		proj.RunSimulation(parameters)
		proj.cSystem.coordinates3 = refcrd3
	#_--------------
	proj.FinishRun() 			
#=====================================================
def MMMD_Heating():
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
	_plotParameters = {"show":True}
	#-------------------------
	parameters = {"protocol":"heating"                  ,
				  "nsteps":5000                         ,
				  "MD_method":"LeapFrog"                ,
				  "temperature_scale_option":"linear"   ,
				  "temperature_scale":  		   15   ,
				  "start_temperature":             20   ,
				  "sampling_factor":               50   ,
				  "log_frequency":				   50   ,
				  "simulation_type":"Molecular_Dynamics",
				  "temperature":               330.15   }
	proj.RunSimulation(parameters)
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
	parameters = {"maxIterations":600 			           ,
				  "rmsGradient":0.1   			           ,
				  "log_frequency":10 		               ,
				  "simulation_type":"Geometry_Optimization",
				  "save_frequency" : 20 		           ,
				  "save_format":".dcd"                     ,
				  "save_pdb": True    			           }
	#---------------------------------------------------
	for alg in algs:
		parameters["optmizer"]=alg
		parameters["trajectory_name"]="QCMMopt_"+alg+".ptGeo"						
		proj.RunSimulation(parameters)
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
	proj=SimulationProject( os.path.join(scratch_path,"QCMM_MDs") )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	#testing qcmm MD 
	#--------------------------------------------------------
	atom1  = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2  = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3  = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")
	atomsf = [ atom1[0], atom2[0], atom3[0] ]
	atom6  = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:O06")
	atom5  = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4  = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")	
	atomss = [ atom4[0], atom5[0], atom6[0] ]
	
	parameters = {"show":True,"ATOMS_RC1":atomsf,"ATOMS_RC2":atomss,"calculate_distances":True,"protocol":"sampling",
				  "nsteps":1000,"sampling_factor":0,"MD_method":"LeapFrog","simulation_type":"Molecular_Dynamics"}
	#--------------------------------------------------------			
	proj.RunSimulation(parameters)
	parameters["nsteps"] = 4000 
	parameters["sampling_factor"] = 100
	proj.RunSimulation(parameters)
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
	atom1  = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2  = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3  = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")
	atomsf = [ atom1[0], atom2[0], atom3[0] ]
	atom6  = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:O06")
	atom5  = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:HE2")
	atom4  = AtomSelection.FromAtomPattern(proj.cSystem,"*:HIE.94:NE2")	
	atomss = [ atom4[0], atom5[0], atom6[0] ]
	#-----------------------------------------------------------------
	parameters = {	"protocol":"sampling"     						,
					"nsteps":500    								,
					"MD_method":"LeapFrog"      					,
				 	"ATOMS_RC1":atomsf          					,
				 	"ATOMS_RC2":atomss          					,
				 	"force_constant_1":100.0    					,
				 	"force_constant_2":95.0     					,
				 	"ndim":2                    					,
				 	"MC_RC1":True               					,
				 	"MC_RC2":True               					,
				 	"simulation_type":"Restricted_Molecular_Dynamics",
				 	"sampling_factor":0       						}
	#------------------------------------------------
	#Run simulation	
	proj.RunSimulation(parameters)
	parameters["sampling_factor"] = 100
	parameters["nstepsns"]        = 2000
	proj.RunSimulation(parameters)
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
	if not name == "Default": _scanFolder = name
	proj=SimulationProject( os.path.join(scratch_path,_scanFolder) )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	#setting atoms for scan
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atomsf = [ atom1[0], atom2[0] ] 
	#--------------------------------------------
	#set parameters for relaxed surface scan
	parameters = { "ATOMS_RC1"  :atomsf		,
				   "dincre_RC1" :_dincre	,
				   "nsteps_RC1" :_nsteps    ,
				   "ndim"       :1    		,
				   "MC_RC1"     :True       ,
				   "save_format":".dcd"     ,
				   "log_frequency":100      ,
				   "simulation_type":"Relaxed_Surface_Scan",
				   "force_constant":4000.0	}	
	proj.RunSimulation(parameters)		
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
	if not name == "Default": _scanFolder = name
	proj=SimulationProject( os.path.join(scratch_path,_scanFolder) )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	#---------------------------------------------
	#setting atoms for scan
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")
	atomsf = [ atom1[0], atom2[0], atom3[0] ] 
	#setting parameters
	parameters = { "ATOMS_RC1":atomsf     ,
				   "dincre_RC1":_dincre   ,
				   "nsteps_RC1":_nsteps   , 
				   "ndim":1               ,
				   "MC_RC1":True          ,
				   "log_frequency":100    ,
				   "save_format":".dcd"   ,
				   "simulation_type":"Relaxed_Surface_Scan",
				   "force_constant":4000.0}
    #run the simulation
    #---------------------------------------------------------------------
	proj.RunSimulation(parameters)		
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
	parameters = { "ATOMS_RC1":atomsf     ,
				   "nsteps_RC1":_nsteps   ,
				   "rc_type_1" :"dihedral", 
				   "ndim":1               ,
				   "MC_RC1":True          ,
				   "save_format":".dcd"   ,
				   "log_frequency":50.0   ,
				   "simulation_type":"Relaxed_Surface_Scan",
				   "force_constant":100.0  }
    #run the simulation
    #---------------------------------------------------------------------
	proj.RunSimulation(parameters)		
	proj.FinishRun()
#=====================================================
def Scan2D_Dihedral(_xnsteps,_ynsteps,name="Default"):
	'''
	Test relaxed scan with dihedral reaction coordinates
	'''
	#---------------------------------------------	
	_scanFolder = "SCAN2D_dihedral"
	if not name == "Default": _scanFolder = name
	proj=SimulationProject( os.path.join(scratch_path,_scanFolder) )		
	proj.LoadSystemFromSavedProject( balapkl )
	proj.cSystem.Summary()
	#---------------------------------------------
	#setting atoms for scan	
	atomsf = [ 4, 6,  8, 14] 
	atomss = [ 6, 8, 14, 16]
	#setting parameters
	parameters = { "ATOMS_RC1":atomsf     ,
				   "ATOMS_RC2":atomss     ,
				   "contour_lines":12     ,
				   "nsteps_RC1":_xnsteps  ,
				   "nsteps_RC2":_ynsteps  ,
				   "rc_type_1" :"dihedral", 
				   "rc_type_2" :"dihedral", 
				   "ndim":2               ,
				   "force_constant_1":40.0,
				   "force_constant_2":40.0,
				   "log_frequency":100    ,
				   "simulation_type":"Relaxed_Surface_Scan",
				   "NmaxThreads":8       }
    #run the simulation
    #---------------------------------------------------------------------
	proj.RunSimulation(parameters)
	proj.SaveProject()		
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
	if not name == "Default": _scanFolder = name
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
	parameters = { "ATOMS_RC1":atomsf	  ,
				   "ATOMS_RC2":atomss	  ,
				   "contour_lines":14     ,
				   "dincre_RC1":_dincrex  ,
				   "dincre_RC2":_dincrey  , 
				   "nsteps_RC1":_xnsteps  ,
				   "nsteps_RC2":_ynsteps  , 
				   "ndim": 2 			  ,
				   "log_frequency":100    ,
				   "simulation_type":"Relaxed_Surface_Scan",				 
				   "NmaxThreads":        4}
	proj.RunSimulation(parameters)		
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
	if not name == "Default": _scanFolder = name
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
				   "contour_lines":15     ,
				   "dincre_RC1":_dincrex  ,
				   "dincre_RC2":_dincrey  , 
				   "nsteps_RC1":_xnsteps  ,
				   "nsteps_RC2":_ynsteps  , 
				   "ndim": 2 			  ,
				   "MC_RC1":True          ,
				   "log_frequency":100    ,
				   "simulation_type":"Relaxed_Surface_Scan",
				   "NmaxThreads":4        }
	#--------------------------------------------------
	proj.RunSimulation(parameters)		
	proj.FinishRun()
#=====================================================
def QCMMScan2DmultipleDistance(_xnsteps,_ynsteps,_dincrex,_dincrey,name="Default"):
	'''
	'''
	#-------------------------------------------------
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMMoptimizations()
	_scanFolder = "QCMM_Scan2D_multiple_distance"
	if not name == "Default": _scanFolder = name
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
	parameters = { "ATOMS_RC1":atomsf	  ,
				   "ATOMS_RC2":atomss	  ,
				   "contour_lines":15     ,
				   "dincre_RC1":_dincrex  ,
				   "dincre_RC2":_dincrey  , 
				   "nsteps_RC1":_xnsteps  ,
				   "nsteps_RC2":_ynsteps  , 
				   "ndim": 2 			  ,
				   "MC_RC1":		True  ,
				   "MC_RC2":		True  ,
				   "log_frequency":100    ,
				   "simulation_type":"Relaxed_Surface_Scan",
				   "NmaxThreads":        4}
				 
	proj.RunSimulation(parameters)		
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
	parameters = { "ATOMS_RC1":atomsf	  ,
				   "ATOMS_RC2":atomss	  ,
				   "dincre_RC1":_dincrex  ,
				   "dincre_RC2":_dincrey  , 
				   "nsteps_RC1":_xnsteps  ,
				   "nsteps_RC2":_ynsteps  , 
				   "ndim": 2 			  ,
				   "contour_lines":15     ,
				   "MC_RC1":		True  ,
				   "MC_RC2":		True  ,
				   "NmaxThreads":        8,
				   "log_frequency":100    ,
				   "simulation_type":"Relaxed_Surface_Scan",
				   "adaptative": True     }
	proj.RunSimulation(parameters)		
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
	if not os.path.exists(_path):
		QCMMScanSimpleDistance(10,0.2,name=_name)
	#-------------------------------------------------
	atom1  = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2  = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atomsf = [ atom1[0], atom2[0] ]	
	rc1 = ReactionCoordinate(atomsf,False,0)	
	rc1.SetInformation(proj.cSystem,0.0)
	#-------------------------------------------------	
	US_parameters = { "ATOMS_RC1":atomsf		       ,
				   "ndim": 1 					       ,
				   "sampling_factor":int(nsteps/10)	   ,
				   "equilibration_nsteps":int(nsteps/2),
				   "production_nsteps":int(nsteps)     ,
				   "source_folder":_path 		       ,
				   "MD_method":"LeapFrog"		       ,
				   "MC_RC1":True				       ,
				   "simulation_type":"Umbrella_Sampling",
				   "NmaxThreads":4 				  }
	#-------------------------------------------------
	#RUN umbrella sampling
	proj.RunSimulation(US_parameters,)	
	#-------------------------------------------------
	#path for the ptRes files
	_path = os.path.join( scratch_path, "FE_simple_distance" )	
	PMF_parameters = { "source_folder":_path  ,
						"xwindows":10         ,
						"ywindows":0          ,
						"crd1_label":rc1.label,
				   		"xnbins":20           ,
				   		"ynbins":0            ,
				   		"simulation_type":"PMF_Analysis",
				   		"temperature":300.15}
	#RUN WHAM, calculate PMF and free energy
	proj.RunSimulation(PMF_parameters)
#=====================================================
def FreeEnergy1DSimpleDistanceOPT(nsteps):	
	'''
	'''
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMM_optimizations()
	proj=SimulationProject( os.path.join(scratch_path,"FE_simple_distance_OPT") )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	#-------------------------------------------------
	_name = "SCAN1D_4FEcalculations_simple_distance"
	_path = os.path.join( os.path.join(scratch_path,_name,"ScanTraj.ptGeo") )
	if not os.path.exists(_path):
		QCMMScanSimpleDistance(10,0.2,name=_name)
	#-------------------------------------------------
	atom1  = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2  = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atomsf = [ atom1[0], atom2[0] ]	
	rc1 = ReactionCoordinate(atomsf,False,0)	
	rc1.SetInformation(proj.cSystem,0.0)
	#-------------------------------------------------
	#RUN umbrella sampling
	US_parameters = { "ATOMS_RC1":atomsf		  ,
				   "ndim": 1 					  ,
				   "sampling_factor":nsteps/10	  ,
				   "equilibration_nsteps":nsteps/2,
				   "production_nsteps":nsteps	  ,
				   "source_folder":_path 		  ,
				   "MD_method":"LeapFrog"		  ,
				   "MC_RC1":True				  ,
				   "optimize":True                ,
				   "simulation_type":"Umbrella_Sampling",
				   "NmaxThreads":4 				  }	
	proj.RunSimulation(US_parameters)	
	#-------------------------------------------------
	#path for the ptRes files
	_path = os.path.join( scratch_path, "FE_simple_distance_OPT" )	
	#RUN WHAM, calculate PMF and free energy
	PMF_parameters = { "source_folder":_path            ,
						"xwindows":10                   ,
						"crd1_label":rc1.label          ,
				   		"xnbins":20                     ,
				   		"simulation_type":"PMF_Analysis",
				   		"temperature":300.15}
	proj.RunSimulation(PMF_parameters)
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
	atomsf = [ atom1[0], atom2[0],atom3[0] ]	
	rc1 = ReactionCoordinate(atomsf,False,0)	
	rc1.SetInformation(proj.cSystem,0.0)	
	_plotParameters = { "contour_lines":15,"xwindows":10,"ywindows":0,"crd1_label":rc1.label}
	#-------------------------------------------------	
	US_parameters = { 'ATOMS_RC1':atomsf		    ,
				   "ndim": 1 						,
				   "sampling_factor":nsteps/10    	,
				   "equilibration_nsteps":nsteps/2 	,
				   "production_nsteps":nsteps		,
				   "source_folder":_path 	    	,
				   "MD_method":"LeapFrog"			,
				   "MC_RC1":True					,
				   "simulation_type":"Umbrella_Sampling",
				   "NmaxThreads":4 				    }
	#-------------------------------------------------
	#Run umbrella sampling 
	proj.RunSimulation(US_parameters)

	_path = os.path.join( scratch_path,"FE_multiple_distance")

	PMF_parameters = { "source_folder":_path,
				   "xnbins":20           ,
				   "xwindows":10         ,
				   "crd1_label":rc1.label,
				   "ynbins":0            ,
				   "simulation_type":"PMF_Analysis",
				   "temperature":300.15  }
	#-------------------------------------------------
	#Run umbrella sampling 
	proj.RunSimulation(PMF_parameters)
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
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")
	atomsf = [ atom1[0], atom2[0],atom3[0] ]	
	rc1 = ReactionCoordinate(atomsf,False,0)	
	rc1.SetInformation(proj.cSystem,0.0)	
	
	_name = "SCAN1D_4US_restart_test"
	QCMMScanMultipleDistance(6,0.2,name=_name)
	_path = os.path.join( scratch_path,_name,"ScanTraj.ptGeo" )
	
	atom1 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.cSystem,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.cSystem,"*:GLU.164:OE2")
	rc1 = ReactionCoordinate(atomsf,False,0)	
	rc1.SetInformation(proj.cSystem,0.0)

	#-------------------------------------------------	
	US_parameters = { "ATOMS_RC1":atomsf		    ,
				   "ndim": 1 						,
				   "sampling_factor":nsteps/10    	,
				   "equilibration_nsteps":nsteps/2 	,
				   "production_nsteps":nsteps		,
				   "source_folder":_path 	    	,
				   "MD_method":"LeapFrog"			,
				   "MC_RC1":True				    ,
				   "restart":True					,
				   "simulation_type":"Umbrella_Sampling",
				   "NmaxThreads":4 				    }
	#-------------------------------------------------
	#Run first umbrella sampling 
	proj.RunSimulation(US_parameters)

	_pathpmf = os.path.join(scratch_path,"UmbrellaSampling_Restart")
	
	_PMFparameters = { "source_folder":_pathpmf,
				   "xnbins":10                 ,
				   "crd1_label":rc1.label      ,
				   "xwindows":6               ,
				   "simulation_type":"PMF_Analysis",
				   "temperature":300.15        }
	#-------------------------------------------------
	#Run umbrella sampling 
	proj.RunSimulation(_PMFparameters)
	input()
	#*************************************************************
	QCMMScanMultipleDistance(10,0.2,name=_name)
	USparameters = { "ATOMS_RC1":atomsf				,
				   "ndim": 1 						,
				   "sampling_factor":nsteps/10    	,
				   "equilibration_nsteps":nsteps/2 	,
				   "production_nsteps":nsteps		,
				   "source_folder":_path 	    	,
				   "MD_method":"LeapFrog"			,
				   "MC_RC1":True					,
				   "restart":True                 ,
				   "simulation_type":"Umbrella_Sampling",
				   "NmaxThreads":8 				    }

	proj.RunSimulation(USparameters)
	_PMFparameters["xwindows"]=10
	proj.RunSimulation(_PMFparameters)
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
	rc1 = ReactionCoordinate(atomsf,False,0)	
	rc1.SetInformation(proj.cSystem,0.0)
	rc2 = ReactionCoordinate(atomss,False,0)	
	rc2.SetInformation(proj.cSystem,0.0)
	#------------------------------------------------------------------
	_name = "SCAN2D_4FEcalculations_simple_distance"
	_path = os.path.join( scratch_path,_name,"ScanTraj.ptGeo")
	if not os.path.exists(_path):
		QCMMScan2DsimpleDistance(6,6,0.2,0.2,name=_name)	

	US_parameters = { "ATOMS_RC1":atomsf		    ,
				   "ATOMS_RC2":atomss				,
				   "ndim": 2 						,
				   "sampling_factor":nsteps/10    	,
				   "equilibration_nsteps":nsteps/2  ,
				   "production_nsteps":nsteps		,
				   "source_folder":_path 			,
				   "MD_method":"LeapFrog"			,
				   "MC_RC1":True					,
				   "MC_RC2":True					,
				   "temperature":300.15             ,
				   "simulation_type":"Umbrella_Sampling",
				   "NmaxThreads":4 					}
	
	proj.RunSimulation(US_parameters)
	_path = os.path.join( scratch_path, "FE_2D_simple_distance")	
	PMFparameters = { "source_folder":_path 		,
				   "xnbins":10           	        ,
				   "xwindows":6                     ,
				   "ywindows":6                     ,
				   "crd1_label":rc1.label           ,
				   "crd2_label":rc2.label           ,
				   "ynbins":10           	        ,
				   "simulation_type":"PMF_Analysis" ,		   
				   "temperature":300.15	            }
	#RUN WHAM, calculate PMF and free energy
	proj.RunSimulation(PMFparameters)
	proj.FinishRun()
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
	rc1 = ReactionCoordinate(atomsf,False,0)	
	rc1.SetInformation(proj.cSystem,0.0)
	rc2 = ReactionCoordinate(atomss,False,0)	
	rc2.SetInformation(proj.cSystem,0.0)

	_name = "SCAN2D_4FEcalculations_mixed_distance"
	_path = os.path.join( scratch_path,_name,"ScanTraj.ptGeo")
	if not os.path.exists(_path):
		QCMMScan2DmixedDistance(6,6,0.2,0.2,name=_name)

	USparameters = { 'ATOMS_RC1':atomsf			,
				   "ATOMS_RC2":atomss			,
				   "ndim": 2 					,
				   "sampling_factor":nsteps/10 	,
				   "equilibration_nsteps":nsteps/2,
				   "production_nsteps":nsteps	,
				   "source_folder":_path 		,
				   "MD_method":"LeapFrog"		,
				   "MC_RC1":True				,
				   "MC_RC2":True				,
				   "restart":True               , 
				   "simulation_type":"Umbrella_Sampling",
				   "NmaxThreads":8 				}
	
	proj.RunSimulation(USparameters)

	_path = os.path.join( scratch_path, "FE_2D_mixed_distance")	
	PMFparameters = { "source_folder":_path ,
					"contour_lines":10   ,
				   "xnbins":10           ,
				   "ynbins":10           ,
				   "xwindows":6          ,
				   "ywindows":6          ,
				   "crd1_label":rc1.label,
				   "crd2_label":rc2.label,
				   "simulation_type":"PMF_Analysis",
				   "temperature":300.15	 }
	#RUN WHAM, calculate PMF and free energy
	proj.RunSimulation(PMFparameters)
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
	
	rc1 = ReactionCoordinate(atomsf,False,0)	
	rc1.SetInformation(proj.cSystem,0.0)
	rc2 = ReactionCoordinate(atomss,False,0)	
	rc2.SetInformation(proj.cSystem,0.0)

	_name = "SCAN2D_4FEcalculations_multiple_distance"
	_path = os.path.join( scratch_path,_name,"ScanTraj.ptGeo")
	if not os.path.exists(_path):
		QCMMScan2DmultipleDistance(6,6,0.2,0.2,name=_name)

	USparameters = { "ATOMS_RC1":atomsf				,
				   "ATOMS_RC2":atomss				,
				   "ndim": 2 						,
				   "sampling_factor":nsteps/10		,
				   "equilibration_nsteps":nsteps/2 	,
				   "production_nsteps":nsteps		,
				   "source_folder":_path 			,
				   "MD_method":"LeapFrog"			,
				   "MC_RC1":True					,
				   "MC_RC2":True					,
				   "simulation_type":"Umbrella_Sampling",
				   "NmaxThreads":4 	     			}
	
	proj.RunSimulation(USparameters)

	_path = os.path.join( scratch_path, "FE_2D_multiple_distance")	
	PMFparameters = { "source_folder":_path ,
				   "xnbins":12           ,
				   "ynbins":12           ,
				   "ywindows":6          ,
				   "xwindows":6          ,
				   "crd1_label":rc1.label,
				   "crd2_label":rc2.label,
				   "simulation_type":"PMF_Analysis",
				   "temperature":300.15	 }
	#RUN WHAM, calculate PMF and free energy
	proj.RunSimulation(PMFparameters)
	proj.FinishRun()
#=====================================================
def FreeEnergyDihedral1D(nsteps):
	'''
	'''	
	proj=SimulationProject( os.path.join(scratch_path,"FE_dihedral_1D") )		
	proj.LoadSystemFromSavedProject( balapkl )	
	#-------------------------------------------------
	_name = "SCAN1D_4FEcalculations_dihedral"
	_path = os.path.join( os.path.join(scratch_path,_name,"ScanTraj.ptGeo") )
	if not os.path.exists(_path ): Scan1D_Dihedral(20, name=_name )
	#-------------------------------------------------	
	atomsf = [4, 6,  8, 14]	
	rc1 = ReactionCoordinate(atomsf,False,0)	
	rc1.SetInformation(proj.cSystem,0.0)

	#-------------------------------------------------	
	parameters = { "ATOMS_RC1":atomsf			  ,
				   "ndim": 1 					  ,
				   "sampling_factor":nsteps/10	  ,
				   "rc_type_1":"dihedral"         ,
				   "equilibration_nsteps":nsteps/2,
				   "force_constant_1":25.0        ,
				   "production_nsteps":nsteps	  ,
				   "source_folder":_path 		  ,
				   "MD_method":"LeapFrog"		  ,
				   "simulation_type":"Umbrella_Sampling",
				   "NmaxThreads":4 				  }
	#-------------------------------------------------
	#RUN umbrella sampling
	proj.RunSimulation(parameters)	
	#-------------------------------------------------
	#path for the ptRes files
	_pathUS = os.path.join( scratch_path, "FE_dihedral_1D" )	
	parameters = { "source_folder":_pathUS			,
				   "xwindows":20                     ,
				   "crd1_label":rc1.label           ,
				   "xnbins":40            			,
				   "ynbins":0                       ,
				   "simulation_type":"PMF_Analysis"	,
				   "temperature":300.15	  			}
	#RUN WHAM, calculate PMF and free energy
	proj.RunSimulation(parameters)
#=====================================================
def FreeEnergyDihedral2D(nsteps):
	'''
	'''
	proj=SimulationProject( os.path.join(scratch_path,"FE_2D_dihedral") )		
	proj.LoadSystemFromSavedProject( balapkl )
	
	atomsf = [ 4, 6,  8, 14 ] 
	atomss = [ 6, 8, 14, 16 ]
	
	rc1 = ReactionCoordinate(atomsf,False,0)	
	rc1.SetInformation(proj.cSystem,0.0)
	rc2 = ReactionCoordinate(atomss,False,0)	
	rc2.SetInformation(proj.cSystem,0.0)
	
	_name = "SCAN2D_4FEcalculations_dihedral"
	_path = os.path.join( scratch_path,_name,"ScanTraj.ptGeo")
	if not os.path.exists(_path): Scan2D_Dihedral(12,12,name=_name)

	parameters = { "ATOMS_RC1":atomsf				,
				   "ATOMS_RC2":atomss				,
				   "ndim": 2 						,
				   "sampling_factor":nsteps/10		,
				   "equilibration_nsteps":nsteps/2 	,
				   "production_nsteps":nsteps		,
				   "source_folder":_path 			,
				   "force_constant_1":25.0          ,
				   "force_constant_2":25.0          ,
				   "rc_type_1":"dihedral"           ,
				   "rc_type_2":"dihedral"           ,
				   "MD_method":"LeapFrog"			,				   
				   "simulation_type":"Umbrella_Sampling",
				   "NmaxThreads":4 	     			}
	
	proj.RunSimulation(parameters)

	_pathUS = os.path.join( scratch_path, "FE_2D_dihedral")	
	parameters = { "source_folder":_pathUS           ,
				   "xnbins":20                       ,
				   "ynbins":20                       ,
				   "simulation_type":"PMF_Analysis"  ,
				   "contour_lines":10                ,
				   "xwindows":12                     ,
				   "ywindows":12                     ,
				   "crd1_label":rc1.label            ,
				   "crd2_label":rc2.label            ,
				   "temperature":300.15	             }
	#RUN WHAM, calculate PMF and free energy
	proj.RunSimulation(parameters)
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
	log_path = os.path.join(scratch_path,"QCMM_SCAN1D_simple_distance.log")
	parameters= {"xsize":30,"type":"1D","log_name":log_path,"crd1_label":rc1_sd.label,"simulation_type":"Energy_Plots"}
	proj.RunSimulation(parameters)
	#===================================================================================
	#2D PES plots 
	if not os.path.exists( os.path.join(scratch_path,"QCMM_Scan2D_simple_distance") ):
		QCMMScan2DsimpleDistance(12,12,0.15,0.15)
	log_path = os.path.join(scratch_path,"QCMM_Scan2D_simple_distance.log")
	parameters= {"xsize":12,"ysize":12,"type":"2D","log_name":log_path,"crd1_label":rc1_sd.label,"crd2_label":rc2_sd.label,"contour_lines":10,"simulation_type":"Energy_Plots"}
	proj.RunSimulation(parameters)
	#------------------------------------------------------------
	if not os.path.exists( os.path.join(scratch_path,"QCMM_Scan2D_mixed_distance") ):
		QCMMScan2DmixedDistance(12,12,0.15,0.15)
	log_path = os.path.join(scratch_path,"QCMM_Scan2D_mixed_distance.log")
	parameters = {"xsize":12,"ysize":12,"type":"2D","log_name":log_path,"simulation_type":"Energy_Plots", "contour_lines":10,"crd1_label":rc1_sd.label,"crd2_label":rc2_sd.label}
	proj.RunSimulation(parameters)
	#------------------------------------------------------------
	if not os.path.exists( os.path.join(scratch_path,"QCMM_Scan2D_multiple_distance") ):
		QCMMScan2DmultipleDistance(12,12,0.15,0.15)
	log_path = os.path.join(scratch_path,"QCMM_Scan2D_multiple_distance.log")
	parameters = {"xsize":12,"ysize":12,"type":"2D","log_name":log_path,"crd1_label":rc1_md.label,"crd2_label":rc2_md.label,"contour_lines":10,"simulation_type":"Energy_Plots"}
	proj.RunSimulation(parameters)
	#===================================================================================
	# 1D Free energy and PMF plots 
	if not os.path.exists( os.path.join(scratch_path,"FE_simple_distance") ):
		FreeEnergy1DSimpleDistance(600)
	log_pathFE  = os.path.join(scratch_path,"FE_simple_distance.log")
	log_pathPMF = os.path.join(scratch_path,"FE_simple_distance.dat")
	parameters      = {"xsize":10,"type":"FE1D","log_name":log_pathFE,"crd1_label":rc1_sd.label,"simulation_type":"Energy_Plots"} # xsize is for windowns size
	proj.RunSimulation(parameters)
	parameters      = {"xsize":20,"type":"WHAM1D","log_name":log_pathPMF,"crd1_label":rc1_sd.label,"simulation_type":"Energy_Plots"} # xsize is for bins size for PMF
	proj.RunSimulation(parameters)
	#------------------------------------------------------------------------------------
	if not os.path.exists( os.path.join(scratch_path,"FE_multiple_distance") ):
		FreeEnergy1DMultipleDistance(600)
	log_pathFE  = os.path.join(scratch_path,"FE_multiple_distance.log")
	log_pathPMF = os.path.join(scratch_path,"FE_multiple_distance.dat")
	parameters      = {"xsize":10,"type":"FE1D","log_name":log_pathFE,"crd1_label":rc1_md.label,"simulation_type":"Energy_Plots"} # xsize is for windowns size
	proj.RunSimulation(parameters)
	parameters      = {"xsize":20,"type":"WHAM1D","log_name":log_pathPMF,"crd1_label":rc1_md.label,"simulation_type":"Energy_Plots"} # xsize is for bins size for PMF
	proj.RunSimulation(parameters)
	#------------------------------------------------------------------------------------
	#====================================================================================
	# 2D Free Energy FES and PMF 
	if not os.path.exists( os.path.join(scratch_path,"FE_2D_simple_distance") ):
		FreeEnergy2DsimpleDistance(1000)
	log_pathFE  = os.path.join(scratch_path,"FE_2D_simple_distance.log")
	log_pathPMF = os.path.join(scratch_path,"FE_2D_simple_distance.dat")
	parameters      = {"xsize":6,"ysize":6,"type":"FE2D","log_name":log_pathFE,"crd1_label":rc1_sd.label,"crd2_label":rc2_sd.label,"contour_lines":10,"simulation_type":"Energy_Plots"} # xsize is for windowns size
	proj.RunSimulation(parameters)
	parameters      = {"xsize":12,"ysize":12,"type":"WHAM2D","log_name":log_pathPMF,"crd1_label":rc1_sd.label,"crd2_label":rc2_sd.label,"contour_lines":10,"simulation_type":"Energy_Plots"} # xsize is for bins size for PMF
	proj.RunSimulation(parameters)
	#-------------------------------------------------------------------------------
	if not os.path.exists( os.path.join(scratch_path,"FE_2D_mixed_distance") ):
		FreeEnergy2DmixedDistance(1000)
	log_pathFE  = os.path.join(scratch_path,"FE_2D_mixed_distance.log")
	log_pathPMF = os.path.join(scratch_path,"FE_2D_mixed_distance.dat")
	parameters      = {"xsize":6,"ysize":6,"type":"FE2D","log_name":log_pathFE,"crd1_label":rc1_sd.label,"crd2_label":rc2_sd.label,"contour_lines":10,"simulation_type":"Energy_Plots"} # xsize is for windowns size
	proj.RunSimulation(parameters)
	parameters      = {"xsize":12,"ysize":12,"type":"WHAM2D","log_name":log_pathPMF,"crd1_label":rc1_sd.label,"crd2_label":rc2_sd.label,"contour_lines":10,"simulation_type":"Energy_Plots"} # xsize is for bins size for PMF
	proj.RunSimulation(parameters)
	#-------------------------------------------------------------------------------
	if not os.path.exists( os.path.join(scratch_path,"FE_2D_multiple_distance") ):
		FreeEnergy2DmultipleDistance(1000)
	log_pathFE = os.path.join(scratch_path,"FE_2D_multiple_distance.log")
	log_pathPMF = os.path.join(scratch_path,"FE_2D_multiple_distance.dat")
	parameters      = {"xsize":6,"ysize":6,"type":"FE2D","log_name":log_pathFE,"crd1_label":rc1_sd.label,"crd2_label":rc2_sd.label,"contour_lines":10,"simulation_type":"Energy_Plots"} # xsize is for windowns size
	proj.RunSimulation(parameters)
	parameters      = {"xsize":12,"ysize":12,"type":"WHAM2D","log_name":log_pathPMF,"crd1_label":rc1_sd.label,"crd2_label":rc2_sd.label,"contour_lines":10,"simulation_type":"Energy_Plots"} # xsize is for bins size for PMF
	proj.RunSimulation(parameters)
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
def ReacCoordSearchers(_type):	
	'''
	Test simulation protocols to determine reaction paths
	'''
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMM_optimizations()
	proj=SimulationProject( os.path.join(scratch_path,"ReactionPathsSearchers") )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	#-------------------------------------------------------------------------------
	#generate initial and final coordinates for NEB 
	#generate trajectory for SAW
	#-------------------------------------------------------------------------------
	_name = "SCAN1D_4NEB_and_SAW"
	_path = os.path.join( os.path.join(scratch_path,_name,"ScanTraj.ptGeo") )
	if not os.path.exists(_path):
		QCMMScanMultipleDistance(16,0.09,name=_name)
	init_path    = os.path.join( _path, "frame0.pkl")
	final_path   = os.path.join( _path, "frame15.pkl")
	saddle_coord = os.path.join( _path, "frame9.pkl") 

	#-------------------------------------------------------------------------------
	if _type == "NEB":
		parameters = {  "init_coord":init_path           					,
						"final_coord":final_path         					,
						"traj_bins":16                   					,
						"refine_methods":["rm1","am1","pm3","pddgpm3","pm6"],
						"RMS_growing_intial_string":1.0  					,
						"simulation_type":"NEB"          					,
						"spring_force_constant":500.0    					,
						"rmsGradient":0.10               					,
						"fixed_terminal_images":True    	                }
	#-------------------------------------------------------------------------------
	elif _type == "NEB_traj":
		parameters = {  "traj_source":_path,
						"traj_bins":16                   					,
						"refine_methods":["rm1","am1","pm3","pddgpm3","pm6"],
						"RMS_growing_intial_string":1.0  					,
						"simulation_type":"NEB"          					,
						"spring_force_constant":500.0    					,
						"rmsGradient":0.3               					,
						"fixed_terminal_images":True                       }
	#------------------------------------------------------------------------------
	elif _type == "SAW":
		parameters = {  "traj_source":_path,
						"traj_bins":16                   					,
						"refine_methods":["rm1","am1","pm3","pddgpm3","pm6"],
						"simulation_type":"SAW"          					,
						"rmsGradient":0.30                                  }
	#------------------------------------------------------------------------------					
	elif _type == "BakerSaddle":
		parameters = {  "saddle_coord":saddle_coord           			    ,
						"simulation_type":"Baker_Saddle"          			,
						"rmsGradient":0.10               					}
	#------------------------------------------------------------------------------
	elif _type == "SteepestDescent":
		parameters = {  "traj_source":_path                                  ,
						"refine_methods":["rm1","am1","pm3","pddgpm3","pm6"] ,
						"simulation_type":"Steep_Path_Searcher"          	 ,
						"rmsGradient":0.10    								 ,
						"function_step":0.025                                ,
						"mass_weighting":True                                ,
						"path_step":2.0          		              		}
	#-------------------------------------------------------------------------------
	proj.RunSimulation(parameters)
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
	
	_name = "SCAN1D_4Refinement"
	_path = os.path.join( os.path.join(scratch_path,_name,"ScanTraj.ptGeo") )
	if not os.path.exists(_path):
		QCMMScanMultipleDistance(30,0.05,name=_name)

	parameters = { "xnbins":30			               ,
				   "ynbins":0			               ,
				   "source_folder":_path                 ,
				   "out_folder":os.path.join(scratch_path, "SMO_EnergyRefinement"),
				   "charge":-3		                    ,
				   "multiplicity":1 	                ,
				   "methods_lists":methods              ,					   
				   "NmaxThreads":8 		                ,
				   "simulation_type":"Energy_Refinement",
				   "crd1_label":rc1_md.label            ,
				   "contour_lines":12                   ,
				   "xlim_list": [-1.2,2.0]              ,
				   "Software":"pDynamo"	}

	proj.RunSimulation(parameters)	
#=====================================================
def pDynamoEnergyRef_2D():
	'''
	Test two dimensinal energy refinement with internal pDynamo quantum chemical methods
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
	#---------------------------------------------
	_name = "SCAN2D_4Refinement"
	_path = os.path.join( os.path.join(scratch_path,_name,"ScanTraj.ptGeo") )
	if not os.path.exists(_path):
		QCMMScan2DmultipleDistance(6,6,0.1,0.1,name=_name)
	#---------------------------------------------
	_plotParameters = {	"show":True              ,
						"crd1_label":rc1_md.label,
						"crd2_label":rc2_md.label,
						"contour_lines":12       ,
						"xlim_list": [-1.2,-0.3] ,
						"ylim_list": [-0.9,-0.2] }
	#---------------------------------------------
	parameters = { "xnbins":6			,
				   "ynbins":6			,
				   "source_folder":_path,
				   "out_folder":os.path.join(scratch_path, "SMO2D_EnergyRefinement"),
				   "charge":-3		    ,
				   "multiplicity":1 	,
				   "methods_lists":methods,					   
				   "NmaxThreads":4		,
				   "crd1_label":rc1_md.label,
				   "crd2_label":rc2_md.label,
				   "contour_lines":12       ,
				   "xlim_list": [-1.2,-0.3] ,
				   "ylim_list": [-0.9,-0.2], 
				   "simulation_type":"Energy_Refinement",
				   "Software":"pDynamo"	}
	#---------------------------------------------
	proj.RunSimulation(parameters)
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
	
	#---------------------------------------------
	parameters = { "xnbins":20			,
				   "ynbins":0			,
				   "mopac_keywords":["grad qmmm","ITRY=5000"] ,
				   "source_folder":_path,
				   "out_folder":os.path.join(scratch_path, "MOPAC_EnergyRefinement"),
				   "charge":-3		    ,
				   "multiplicity":1 	,
				   "methods_lists":methods,					   
				   "NmaxThreads":1 		,
				   "crd1_label":rc1_md.label,
				   "contour_lines":12       ,
				   "xlim_list": [-1.2,2.0] ,
				   "simulation_type":"Energy_Refinement",
				   "Software":"mopac"	}
	#---------------------------------------------
	proj.RunSimulation(parameters)	
#=====================================================
def Change_QC_Region():
	pass
#=====================================================
def CombinedFES_ABinitioSMO():
	pass
#=====================================================
def ORCAEnergy():
	pass
#=====================================================
def Thermodynamics():
	pass
#=====================================================
if __name__ == "__main__":
	#MMMD_Algorithms()                         			#TESTED
	#MMMD_Heating()										#TESTED
	#QCMM_Energies()									#TESTED
	#QCMM_DFTBplus()									#TESTED
	#QCMM_Orca()										#TESTED
	#QCMM_optimizations()								#TESTED
	#QCMM_MD()											#TESTED
	#QCMM_MDrestricted()								#TESTED
	#QCMMScanSimpleDistance(20,0.06)					#TESTED
	#QCMMScanMultipleDistance(20,0.06)					#TESTED
	#QCMMScan2DsimpleDistance(10,10,0.2,0.2)			#TESTED
	#QCMMScan2DmixedDistance(10,10,0.2,0.2)				#TESTED
	#QCMMScan2DmultipleDistance(10,10,0.2,0.2)			#TESTED
	#QCMMScans2D_Adaptative(10,10,0.2,0.2)				#TESTED
	#Scan1D_Dihedral(36)								#TESTED
	#Scan2D_Dihedral(10,10)							 	#TESTED
	#FreeEnergy1DSimpleDistance(500)
	#FreeEnergy1DMultipleDistance(500)
	#FreeEnergyDihedral1D(2000)
	#FreeEnergy1DSimpleDistanceOPT(500)
	#UmbrellaSampling1Drestart(500)
	#FreeEnergy2DsimpleDistance(500)	
	#FreeEnergy2DmixedDistance(500)		
	#FreeEnergy2DmultipleDistance(500)
	#pDynamoEnergyRef_1D()								#TESTED
	#EnergyAnalysisPlots()								#TESTED
	ReacCoordSearchers("BakerSaddle")					#TESTED
	#MopacEnergyRef()									#TESTED
	#pDynamoEnergyRef_2D()								#TESTED

