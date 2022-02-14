#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = SimulationsPreset.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#--------------------------------------------------------------
import os, glob
#--------------------------------------------------------------
#Loading own libraries
from commonFunctions 		import *
from LogFile  				import LogFile
from GeometrySearcher 		import GeometrySearcher
from RelaxedScan 			import SCAN
from MolecularDynamics  	import MD
from UmbrellaSampling  		import US
from PotentialOfMeanForce 	import PMF
from ReactionCoordinate 	import *
from EnergyRefinement	 	import *
from Analysis    			import EnergyAnalysis
from Analysis 				import DistanceAnalysis

#--------------------------------------------------------------
#loading pDynamo Libraries
from pBabel                    import *                                     
from pCore                     import *
#---------------------------------------                                     
from pMolecule                 import *                              
from pMolecule.MMModel         import *
from pMolecule.NBModel         import *                                     
from pMolecule.QCModel         import *
#---------------------------------------
from pScientific               import *                                     
from pScientific.Arrays        import *                                     
from pScientific.Geometry3     import *                                     
from pScientific.RandomNumbers import *                                     
from pScientific.Statistics    import *
from pScientific.Symmetry      import *
#--------------------------------------                                    
from pSimulation               import *
#=============================================================
class Simulation:
	'''
	Class to set up preset simulations to be perfomed
	'''
	def __init__(self,_system,_simulationType,baseFolder):
		'''
		'''
		self.molecule 			= _system
		self.simulationType 	= _simulationType
		self.baseFolder  		= baseFolder# the baseFolder for the simulations will be the current dir for now.
		self.MAXnprocs 			= 1 # maximum number of virtual threads to be used in the simulations
		self.coorddinatesFolder = "" # Name of the folder containing the pkls to be read. Used in more than one preset here
		self.logFreq 			= 1
		self.samplingFactor 	= 1 # this is usually let to the default class value, unless the user want to modify
		self.nProcs 			= NmaxThreads
		#for restricted simulations
		
		#enviromental parameters and their default values 
		self.temperature = 300.15
		self.pressure 	 = 1

		#specific parameters for molecular dynamics run
		self.equiNsteps = 0
		self.prodNsteps = 0	

		#--------------------------------------------------
		if not os.path.exists( self.baseFolder ):
			os.makedirs(self.baseFolder)

	#=======================================================================
	def Execute(self,_parameters,_plotParameters=None):
		'''
		Function to call the class method to execute the preset simulation
		Parameters:
			_parameters    : python dict with parameters for simulation
			_plotParameters: python dict with parameters for plot graphics and analysis
		'''		
		#-------------------------------------------------------------
		if self.simulationType == "Energy_Refinement":			
			self.EnergyRefine(_parameters,_plotParameters)		
		#-------------------------------------------------------------
		elif self.simulationType == "Geometry_Optimization":			
			self.GeometryOptimization(_parameters)
		#-------------------------------------------------------------
		elif self.simulationType == "Relaxed_Surface_Scan":			
			self.RelaxedSurfaceScan(_parameters,_plotParameters)
		#-------------------------------------------------------------
		elif self.simulationType == "Molecular_Dynamics":
			self.MolecularDynamics(_parameters,_plotParameters)
		#-------------------------------------------------------------	
		elif self.simulationType == "Restricted_Molecular_Dynamics":			
			self.RestrictedMolecularDynamics(_parameters,_plotParameters)
		#-------------------------------------------------------------
		elif self.simulationType == "Umbrella_Sampling":
			self.UmbrellaSampling(_parameters,_plotParameters)
		#-------------------------------------------------------------
		elif self.simulationType == "PMF_Analysis":
			self.PMFAnalysis(_parameters,_plotParameters)		
		#-------------------------------------------------------------
		elif self.simulationType == "Normal_Modes":				
			self.NormalModes(_parameters,_plotParameters)
		#-------------------------------------------------------------
		elif self.simulationType == "Delta_Free_Energy":			
			self.DeltaFreeEnergy(_parameters)
		#-------------------------------------------------------------
		elif self.simulationType == "NEB":
			self.NEB(_parameters)
		#-------------------------------------------------------------
		elif self.simulationType == "SAW":
			self.SAW(_parameters)
		#-------------------------------------------------------------
		elif self.simulationType == "Simulating_Annealing":
			self.SimulatingAnnealing(_parameters)
		#-------------------------------------------------------------
		elif self.simulationType == "Steered_Molecular_Dynamics":
			self.SMD(_parameters)
		#-------------------------------------------------------------
		elif self.simulationType == "Trajectory_Analysis":
			self.TrajectoryAnalysis(_parameters,_plotParameters) 
		#-------------------------------------------------------------
		elif self.simulationType == "Energy_Plots":
			self.EnergyPlots(_parameters,_plotParameters)

	#==================================================================
	def EnergyRefine(self,_parameters,_plotParameters):
		'''
		Set up and execute energy refinement using a series of methods
		Parameters:
			_parameters: python dict with parameters for simulation
		'''
		_Restart      = False
		dimensions    = [ 0,0 ] 
		dimensions[0] =  _parameters["xnbins"] 
		if "ynbins" in _parameters:
			dimensions[1] = _parameters["ynbins"]
		if "restart" in _parameters:
			_Restart = True
		#------------------------------------------------------------------
		ER = EnergyRefinement(self.molecule  					,
							  _parameters["source_folder"]  	,
							  _parameters["out_folder"]         ,dimensions,
							  _parameters["charge"]             ,
							  _parameters["multiplicity"]		)
		#------------------------------------------------------------------
		if "change_qc_region" in _parameters:
			ER.ChangeQCRegion(_parameters["center"],_parameters["radius"])
			#------------------------------------------------------------
		if _parameters["Software"] == "pDynamo":
			ER.RunInternalSMO(_parameters["methods_lists"],_parameters["NmaxThreads"])
			#------------------------------------------------------------
		elif _parameters["Software"] == "DFTBplus":
			pass
			#------------------------------------------------------------
		elif _parameters["Software"] == "Mopac":
			pass
			#------------------------------------------------------------
		elif _parameters["Software"] == "ORCA":
			ER.RunORCA(_parameters["orca_method"],_parameters["basis"],_parameters["NmaxThreads"],_restart=_Restart)
			#------------------------------------------------------------
				
		#===========================================================
		#Set plor parameters
		cnt_lines  = 12
		crd1_label = "Reaction Coordinate #1"
		crd2_label = "Reaction Coordinate #2"
		xlims = [ 0, dimensions[0] ]
		ylims = [ 0, dimensions[1] ]
		show  = False
		#check parameters for plot
		if "contour_lines" in _plotParameters:
			cnt_lines  = _plotParameters["contour_lines"]
		if "crd1_label" in _plotParameters:
			crd1_label = _plotParameters["crd1_label"]
		if "crd2_label" in _plotParameters:
			crd2_label = _plotParameters["crd2_label"]
		if "xlim_list" in _plotParameters:
			xlim = _plotParameters["xlim_list"]
		if "ylim_list" in _plotParameters:
			ylim = _plotParameters["ylim_list"]
		if "show" in _plotParameters:
			show = True
		#------------------------------------------------------------
		ER.WriteLog()
		if len(dimensions[1]) > 0:
			TYPE = "2DRef"
		else: 
			TYPE = "1DRef"		
		EA = EnergyAnalysis(dimensions[0],dimensions[1],_type=TYPE)
		EA.ReadLog( os.path.join(ER.baseName,"EnergyRefinement.log") )
		#-------------------------------------------------------------
		if len(dimensions[1]) > 0:
			EA.Plot2D(cnt_lines,crd1_label,crd2_label,xlim,ylim,show)
		else:
			if "methods_lists" in _parameters:
				if len(_parameters["methods_lists"]) > 1:
					EA.MultPlot1D(_plotParameters["crd1_label"])
			else:
				EA.Plot1D(_plotParameters["crd1_label"],show)
	#==================================================================
	def GeometryOptimization(self,_parameters):
		'''
		Set up and execture the search of local minima for the system passed
		Parameters:
			_parameters: python dict with parameters for simulation
		'''
		_Optimizer = "ConjugatedGradient"
		if "optmizer" in _parameters:
			_Optimizer = _parameters["optmizer"]

		Gopt = GeometrySearcher(self.molecule,self.baseFolder)		
		Gopt.ChangeDefaultParameters(_parameters)
		Gopt.Minimization(_Optimizer)
		Gopt.Finalize()
	#==================================================================
	def RelaxedSurfaceScan(self,_parameters,_plotParameters):
		'''
		Set up and execute one/two-dimensional relaxed surface scans 
		Parameters:
			_parameters: python dict with parameters for simulation
			_plotParameters: python dict with parameters for plot graphics and analysis
		'''
		#------------------------------------------------------------------
		_Adaptative = False
		_Optmizer   = "ConjugatedGradient"
		MCR1 		= False
		MCR2 		= False
		RD          = _parameters["ndim"]
		rcType1     = "Distance"
		rcType2     = "Distance"
		nDims       = _parameters['ndim']
		#-------------------------------------------------------------------
		if "optmizer" in _parameters:
			_Optmizer   = _parameters["optmizer"]
		if "adaptative" in _parameters:
			_Adaptative = True
		if "MC_RC1" in _parameters:
			MCR1 = True
		if "MC_RC2" in _parameters:
			MCR2 = True	
		#--------------------------------------------------------------------
		scan = SCAN(self.molecule,self.baseFolder,_Optmizer,ADAPTATIVE=_Adaptative)
		scan.ChangeDefaultParameters(_parameters)	
		#--------------------------------------------------------------------
		rc1 = ReactionCoordinate(_parameters["ATOMS_RC1"], MCR1,_type=rcType1)
		rc1.SetInformation(self.molecule,_parameters['dincre_RC1'])
		rc2 = None
		if nDims == 2:
			rc2 = ReactionCoordinate(_parameters["ATOMS_RC2"], MCR2,_type=rcType2)
			rc2.SetInformation(self.molecule,_parameters['dincre_RC2'])				
		#------------------------------------------------------
		scan.SetReactionCoord(rc1)
		xlims = [0, _parameters['nSteps_RC1'] ]
		ylims = [0,0]
		#--------------------------------------------------------------------------------
		if nDims == 2:
			scan.SetReactionCoord(rc2)
			scan.Run2DScan(_parameters['nSteps_RC1'], _parameters['nSteps_RC2'] )
			ylims = [ 0,  _parameters['nSteps_RC2']]
		elif nDims == 1:
			scan.RunONEDimensionSCAN(_parameters['nSteps_RC1'])
		#...............
		scan.Finalize()		
		#================================================================
		#Set plor parameters
		cnt_lines  = 12
		crd1_label = rc1.label
		crd2_label = ""
		nRC2 = 0
		show = False
		if nDims == 2:
			crd2_label = rc2.label
		xlims = [ 0,  _parameters['nSteps_RC1'] ]
		#check parameters for plot
		if "contour_lines" in _plotParameters:
			cnt_lines  = _plotParameters["contour_lines"]		
		if "xlim_list" in _plotParameters:
			xlims = _plotParameters["xlim_list"]
		if "ylim_list" in _plotParameters:
			ylims = _plotParameters["ylim_list"]
		if "nSteps_RC2" in _parameters:
			nRC2= _parameters["nSteps_RC2"]
		if "show" in _plotParameters:
			show = True
		#------------------------------------------------------------
		if nDims == 2:
			TYPE = "2D"
		elif nDims == 1: 
			TYPE = "1D"		

		#------------------------------------------------------------
		EA = EnergyAnalysis(_parameters['nSteps_RC1'],nRC2,_type=TYPE)
		EA.ReadLog( scan.baseName+"_SCAN{}D.log".format(nDims)) 
		#-------------------------------------------------------------
		if nDims == 2:
			EA.Plot2D(cnt_lines,crd1_label,crd2_label,xlims,ylims,show)
		elif nDims == 1:
			EA.Plot1D(crd1_label,show)
		#-------------------------------------------------------------
	#=================================================================
	def MolecularDynamics(self,_parameters,_plotParameters):
		'''
		Set up and execute molecular dynamics simulations.
		Parameters:
			_parameters: python dict with parameters for simulation
			_plotParameters: python dict with parameters for plot graphics and analysis
		'''
		#-------------------------------------------------------------
		MDrun = MD(self.molecule,self.baseFolder,_parameters['MD_method'])
		
		#If there is some key in _parameters set to modify some varibles then this is done here
		#If there is any this next command will have no effect
		MDrun.ChangeDefaultParameters(_parameters)
		#---------------------------------------------------------------
		if "protocol" in _parameters:
			if _parameters["protocol"] 		== "heating":
				print(_parameters['production_nsteps'])
				MDrun.HeatingSystem(_parameters['production_nsteps'])
			elif _parameters["protocol"] 	== "equilibration":
				MDrun.RunEquilibration(_parameters['equilibration_nsteps'])
			elif _parameters["protocol"] 	== "production":
				MDrun.RunEquilibration(_parameters['equilibration_nsteps'])
				MDrun.RunProduction(_parameters['production_nsteps'])
		#----------------------------------------------------------------
		if not _plotParameters == None:
			show = False
			RCs  = None
			if "show" in _plotParameters:
				show = True
			t_time = _parameters["production_nsteps"]*0.001
			DA = DistanceAnalysis(MDrun.trajectoryNameCurr,self.molecule,t_time)
			DA.CalculateRG_RMSD()
			DA.PlotRG_RMS(show)
			
			if "calculate_distances" in _plotParameters:
				rc1 = ReactionCoordinate(_plotParameters["ATOMS_RC1"],False,0)
				rc1.SetInformation(self.molecule,0)
				RCs = [rc1]
				rc2 = None
				if "ATOMS_RC2" in _plotParameters:
					rc2 = ReactionCoordinate(_plotParameters["ATOMS_RC2"],False,0)
					rc2.SetInformation(self.molecule,0)
					RCs.append(rc2)
				DA.DistancePlots()
	#==================================================================
	def RestrictedMolecularDynamics(self,_parameters,_plotParameters):
		'''
		Set up and execute molecular dynamics simulations.
		Parameters:
			_parameters: python dict with parameters for simulation
			_plotParameters: python dict with parameters for plot graphics and analysis
		'''
		#----------------------------------------------------------------
		restraints = RestraintModel( )
		self.molecule.DefineRestraintModel( restraints )
		
		MCR1 = False
		MCR2 = False
		if "MC_RC1" in _parameters:
			MCR1 = True
		if "MC_RC2" in _parameters:
			MCR2 = True
		#---------------
		rcType1 = "Distance"
		rcType2 = "Distance"
		if "type_rc1" in _parameters:
			rcType1 = _parameters["type_rc1"]
		if "type_rc2" in _parameters:
			rcType2 = _parameters["type_rc2"]

		#-------------------------------------------------------------------
		forcK = _parameters["forceC"]		
		restrainDimensions = _parameters['ndim']
		#-------------------------------------------------------------------
		rc1 = ReactionCoordinate(_parameters["atoms_M1"],MCR1,_type=rcType1)
		rc1.SetInformation(self.molecule,0)
		nDims = _parameters['ndim']
		rc2 = None
		if nDims == 2:
			rc2 = ReactionCoordinate(_parameters["atoms_M2"],MCR2,_type=rcType2)
			rc2.SetInformation(self.molecule,0)
		#-------------------------------------------------------------------
		distance = rc1.minimumD
		rmodel = RestraintEnergyModel.Harmonic( distance, forcK )
		if rc1.nAtoms == 3:				
			restraint = RestraintMultipleDistance.WithOptions( energyModel=rmodel, distances=[ [ rc1.atoms[1], rc1.atoms[0], rc1.weight13 ], [ rc1.atoms[1], rc1.atoms[2], rc1.weight31 ] ] ) 
		elif rc1.nAtoms == 2:				
			restraint = RestraintDistance.WithOptions( energyModel=rmodel, point1=rc1.atoms[0], point2=rc1.atoms[1] )
		restraints['M1'] =  restraint
		#-------------------------------------------------------------------
		if nDims == 2:
			distance = rc2.minimumD
			rmodel = RestraintEnergyModel.Harmonic( distance, forcK )
			if rc2.nAtoms == 3:				
				restraint = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ rc2.atoms[1], rc2.atoms[0], rc2.weight13 ], [ rc2.atoms[1], rc2.atoms[2], rc2.weight31 ] ] ) 
			elif rc1.nAtoms == 2:				
				restraint = RestraintDistance.WithOptions( energyModel=rmodel, point1=rc2.atoms[0], point2=rc2.atoms[1] )
			restraints['M2'] =  restraint		
		#----------------------------------------------------------------
		MDrun = MD(self.molecule,self.baseFolder,_parameters['MD_method'])
		MDrun.ChangeDefaultParameters(_parameters)
		MDrun.RunProductionRestricted(_parameters['equilibration_nsteps'],_parameters['production_nsteps'],_parameters["sampling_factor"])
		#-----------------------------------------------------------------		
		if not _plotParameters == None:
			t_time = _parameters["production_nsteps"]*0.001
			if "show" in _plotParameters:
				show = True
			DA = DistanceAnalysis(MDrun.trajectoryNameCurr,self.molecule,t_time)
			DA.CalculateRG_RMSD()
			DA.PlotRG_RMS(show)				
			RCs = [rc1]
			if nDims > 1:				
				RCs.append(rc2)							
			DA.DistancePlots(RCs,show)
			DA.ExtractFrames()
	#=======================================================================
	def UmbrellaSampling(self,_parameters,_plotParameters):
		'''
		Set up and execute umbrella sampling simulations and Free energy calculations for reaction path trajectory.
		Parameters:
			_parameters: python dict with parameters for simulation
			_plotParameters: python dict with parameters for plot graphics and analysis
		'''
		#---------------------------------------
		MCR1 = False
		MCR2 = False		
		rcType1 = "Distance"
		rcType2 = "Distance"
		#---------------------------------------
		if "MC_RC1" in _parameters:
			MCR1 = True
		if "MC_RC2" in _parameters:
			MCR2 = True
		#---------------------------------------
		_Restart 	= False
		_Adaptative = False
		if "restart" in _parameters:
			_Restart = True 
		if "adaptative" in _parameters:
			_Adaptative = True
		#-------------------------------------------------------------------
		rc1 = ReactionCoordinate(_parameters["atoms_M1"],MCR1,_type=rcType1)
		rc1.SetInformation(self.molecule,0)
		nDims = _parameters['ndim']
		rc2 = None
		if nDims == 2:
			rc2 = ReactionCoordinate(_parameters["atoms_M2"],MCR2,_type=rcType2)
			rc2.SetInformation(self.molecule,0)
		#---------------------------------------
		USrun = US(self.molecule  						,
			       self.baseFolder 						,
			       _parameters['equilibration_nsteps']  ,
			       _parameters['production_nsteps']     ,
			       _parameters["MD_method"]             ,
			       RESTART=_Restart                     ,
			       ADAPTATIVE=_Adaptative               )
		#---------------------------------------
		USrun.ChangeDefaultParameters(_parameters)
		USrun.SetMode(rc1)
		if _parameters["ndim"] == 1:
			USrun.Run1DSampling(_parameters["source_folder"],_parameters["sampling_factor"])
		elif _parameters["ndim"] == 2:
			USrun.SetMode(rc2)
			USrun.Run2DSampling(_parameters["source_folder"],_parameters["sampling_factor"])		
	#=========================================================================
	def PMFAnalysis(self,_parameters,_plotParameters):
		'''
		Calculate potential of mean force and Free energy from restricted molecular dynamics
		Parameters:
			_parameters: python dict with parameters for simulation
			_plotParameters: python dict with parameters for plot graphics and analysis
		'''
		potmean = PMF( self.molecule, _parameters["source_folder"], self.baseFolder )
		potmean.CalculateWHAM(_parameters["xnbins"],_parameters["ynbins"],_parameters["temperature"])
	#=========================================================================
	def NormalModes(self,_parameters):
		'''
		Simulation preset to calculate the normal modes and to write thr trajectory for a specific mode.
		Parameters:
			_parameters: python dict with parameters for simulation
		'''	
		mode 		= 0
		temperature = 300.15
		Cycles 		= 10 
		Frames  	= 10 
		if "temperature" in _parameters:
			temperature = _parameters['temperature']
		if "cycles" in _parameters:
			Cycles = _parameters['cycles']
		if "frames" in _parameters:
			Frames = _parameters['frames']
		if "mode" in _parameters:
			mode = _parameters["mode"]


		NormalModes_SystemGeometry ( self.molecule, modify = ModifyOption.Project )
		if _mode > 0:
			trajectory = ExportTrajectory ( os.path.join (self.baseFolder, "NormalModes","trj"), self.molecule )
			NormalModesTrajectory_SystemGeometry(	self.molecule		      ,
                                       			 	trajectory                ,
                                       				mode        = _mode	      ,
                                       				cycles      = Cycles      ,
                                       				frames      = Frames 	  ,
                                       				temperature = temperature )
	#==========================================================================
	def DeltaFreeEnergy(self,_parameters):
		'''
		Calculate the free energy difference between two configurations of the system using the 
		statistical thermodynamics partition functions from through the normal modes calculations
		Parameters:
			_parameters: python dict with parameters for simulation
		'''
		
		#initial Structure
		pressure       = 1.0
		temperature    = 300.15 
		symmetryNumber = 1

		if "pressure" in _parameters:
			pressure = _parameters["pressure"]

		self.molecule.coordinates3 = ImportCoordinates3(_parameters["initial_coordinates"])
		e0 = self.molecule.Energy()
		NormalModes_SystemGeometry( self.molecule, modify = ModifyOption.Project )
		Gibbs = [] 
		tdics = ThermodynamicsRRHO_SystemGeometry ( self.molecule 							,
                                                    pressure       = pressure       	,
                                                    symmetryNumber = self.symmetryNumber 	,
                                                    temperature    = self.temperature    	)
		Gibbs.append( tdics["Gibbs Free Energy"] )
    	# Final struct
		self.molecule.coordinates3 = ImportCoordinates3(_parameters["final_coordinates"])
		e1 = self.molecule.Energy()
		NormalModes_SystemGeometry ( self.molecule, modify = ModifyOption.Project )
    	 
		tdics = ThermodynamicsRRHO_SystemGeometry ( self.molecule 							,
                                                    pressure       = self.pressure       	,
                                                    symmetryNumber = self.symmetryNumber 	,
                                                    temperature    = self.temperature    	)
		Gibbs.append( tdics["Gibbs Free Energy"] )
	#=========================================================================	
	def NEB(self,_parameters,_plotParameters):
		'''
		Class method to set up and execute Nudget Elastic Band simulations to generate a reaction path trajectory
		Parameters:
			_parameters: python dict with parameters for simulation
			_plotParameters: python dict with parameters for plot graphics and analysis
		'''
		NEBrun = GeometrySearcher(self.system,self.baseFolder)

		#if there any parameters to be modified 
		NEBrun.ChengeDefaultParameters(_parameters)
		NEBrun.NudgedElasticBand(_parameters['initial_coordinates'], 	\
								 _parameters['final_coordinates'] ,  	\
								 _parameters['NEB_nbins'],				\
								 _parameters["RMS_growing_intial_string"] )
	#=========================================================================
	def SAW(self,_parameters):
		'''
		Set up and execute Self-Avoid Walking simulations to generate a reaction path trajectory 
		'''
	#=========================================================================
	def SimulatingAnnealing(self,_parameters):
		'''
		Set up and execute Simulate annealing simulations		
		'''
		pass
	#=========================================================================
	def SMD(self,_parameters):
		'''
		Set up and execute Steered Molecular Dynamics simulations
		'''
		pass
#=============================================================================
#========================END OF THE FILE======================================
#=============================================================================