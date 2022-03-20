 #!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = SimulationsPreset.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#--------------------------------------------------------------
import os, glob, sys
#--------------------------------------------------------------
VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')
if not VISMOL_HOME == None:
	sys.path.append(os.path.join(VISMOL_HOME,"easyhybrid/pDynamoMethods") ) 
else:
	sys.path.append(os.path.join("/home/igorchem/easyhybrid/pDynamoMethods") ) 
#Loading own libraries
#-------------------------------------------------------------
from EnergyAnalysis     	import EnergyAnalysis
from TrajectoryAnalysis 	import TrajectoryAnalysis
#-------------------------------------------------------------
from GeometrySearcher 	    import GeometrySearcher
from RelaxedScan 			import SCAN
from MolecularDynamics  	import MD
from UmbrellaSampling  	    import US
from PotentialOfMeanForce   import PMF
from ReactionCoordinate 	import ReactionCoordinate
from EnergyRefinement	 	import EnergyRefinement
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
		Deafault constructor
		'''
		self.molecule 		= _system
		self.simulationType = _simulationType
		self.baseFolder  	= baseFolder# the baseFolder for the simulations will be the current dir for now.
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
		#-------------------------------------------------------------------------------------------------------------------------
		if 	 self.simulationType == "Energy_Refinement": 				self.EnergyRefine(_parameters,_plotParameters)		
		elif self.simulationType == "Geometry_Optimization":			self.GeometryOptimization(_parameters)
		elif self.simulationType == "Relaxed_Surface_Scan":	 			self.RelaxedSurfaceScan(_parameters,_plotParameters)
		elif self.simulationType == "Molecular_Dynamics":				self.MolecularDynamics(_parameters,_plotParameters)	
		elif self.simulationType == "Restricted_Molecular_Dynamics": 	self.RestrictedMolecularDynamics(_parameters,_plotParameters)
		elif self.simulationType == "Umbrella_Sampling":				self.UmbrellaSampling(_parameters,_plotParameters)
		elif self.simulationType == "PMF_Analysis":						self.PMFAnalysis(_parameters,_plotParameters)			
		elif self.simulationType == "Normal_Modes":						self.NormalModes(_parameters,_plotParameters)		
		elif self.simulationType == "Delta_Free_Energy":				self.DeltaFreeEnergy(_parameters)		
		elif self.simulationType == "NEB":								self.NEB(_parameters,_plotParameters)		
		elif self.simulationType == "SAW":								self.SAW(_parameters,_plotParameters)		
		elif self.simulationType == "Simulating_Annealing":				self.SimulatingAnnealing(_parameters)		
		elif self.simulationType == "Steered_Molecular_Dynamics":		self.SMD(_parameters)		
		elif self.simulationType == "Trajectory_Analysis":				self.TrajectoryPlots(_parameters,_plotParameters) 		
		elif self.simulationType == "Energy_Plots":						self.EnergyPlots(_parameters,_plotParameters)
	#=================================================================================================================
	def EnergyRefine(self,_parameters,_plotParameters):
		'''
		Set up and execute energy refinement using a series of methods
		Parameters:
			_parameters: python dict with parameters for simulation
				Mandatory:
					"xbins"			: Number of frames for first/only coordinate 
					"source_folder" : path of folder containing frames to refine 
					"out_folder"    : path to output logs and other results
					"charge"        : charge for QM region
					"multiplicity"  : multiplicity for QM region
					"Software"  	: engine used to calculate the energy refinement
				Optinal :
					"ybins" 		  : Number of frames for second coordinate
					"change_qc_region": Flag indicating the intention of modifying the QC regions
					"center"		  : The center of the new QC region
					"radius"		  : The radius from the center of the new QC region
					"orca_method"     : Energy method provided in ORCA (eg.: HF, b3lyp, mp2...)
					"basis"           : String containing the orca key for basis functions
					"NmaxThreads"     : Maximum number of virtual threads to be used by pymp
			_plotParameters:python dict with paramters for post-analysis and plots
				Mandatory :
				Optinal   : 
		'''
		_Restart      = False
		dimensions    = [0,0] 
		dimensions[0] =  _parameters["xnbins"]
		nmaxthreads   = 1 
		if "ynbins"  in _parameters:  dimensions[1] = _parameters["ynbins"]
		if "restart" in _parameters: _Restart = True
		if "NmaxThreads" in _parameters: nmaxthreads = _parameters["NmaxThreads"]
		#------------------------------------------------------------------
		ER = EnergyRefinement(self.molecule  					,
							  _parameters["source_folder"]  	,
							  _parameters["out_folder"]         ,dimensions,
							  _parameters["charge"]             ,
							  _parameters["multiplicity"]		)
		#------------------------------------------------------------------
		if "change_qc_region" in _parameters:	 	ER.ChangeQCRegion(_parameters["center"],_parameters["radius"])
		if 	 _parameters["Software"] == "pDynamo":  ER.RunInternalSMO(_parameters["methods_lists"],nmaxthreads)
		elif _parameters["Software"] == "DFTBplus":	ER.RunDFTB()
		elif _parameters["Software"] == "mopac" or _parameters["Software"]=="MOPAC":
			_mopacKeyWords = ["AUX","LARGE"] 
			if "mopac_keywords" in _parameters:
				for key in _parameters["mopac_keywords"]:
					_mopacKeyWords.append(key)
			ER.RunMopacSMO(_parameters["methods_lists"],_mopacKeyWords)
			#------------------------------------------------------------
		elif _parameters["Software"] == "ORCA":
			ER.RunORCA(_parameters["orca_method"],_parameters["basis"],nmaxthreads,_restart=_Restart)
		#===============================================================
		#Set plor parameters
		cnt_lines  	= 12
		crd1_label 	= "Reaction Coordinate #1"
		crd2_label 	= "Reaction Coordinate #2"
		xlim 		= [ 0, dimensions[0] ]
		ylim 		= [ 0, dimensions[1] ]
		show  		= False
		#check parameters for plot
		if "contour_lines" in _plotParameters: cnt_lines  = _plotParameters["contour_lines"]
		if "crd1_label"    in _plotParameters: crd1_label = _plotParameters["crd1_label"]
		if "crd2_label"    in _plotParameters: crd2_label = _plotParameters["crd2_label"]
		if "xlim_list"     in _plotParameters: xlim       = _plotParameters["xlim_list"]
		if "ylim_list"     in _plotParameters: ylim 	  = _plotParameters["ylim_list"]
		if "show" 		   in _plotParameters: show       = _plotParameters["show"]
		#------------------------------------------------------------
		ER.WriteLog()
		if dimensions[1] > 0: TYPE = "2DRef"
		else: TYPE = "1DRef"		
		EA = EnergyAnalysis(dimensions[0],dimensions[1],_type=TYPE)
		EA.ReadLog( os.path.join(ER.baseName+"_energy.log") )
		#-------------------------------------------------------------
		if dimensions[1] > 0:
			EA.MultPlot2D(cnt_lines,crd1_label,crd2_label,xlim,ylim,show)
		else:
			if "methods_lists" in _parameters:
				if len(_parameters["methods_lists"]) > 1: EA.MultPlot1D(crd1_label)
			else: EA.Plot1D(crd1_label,show)
	#==================================================================
	def GeometryOptimization(self,_parameters):
		'''
		Set up and execture the search of local minima for the system passed
		Parameters:		
			_parameters: python dict with parameters for simulation
				Mandatory ( if not provided a "key-Error" will be thrown ):
					optimizer: name of the optimization algorithm					 
				Optinal :
					trajectory_name: name to save the trajectory
					maxIterations: maximum number of itetarions (integer) 
					log_frequency: log frequency  (integer)
					save_pdb     : whether to save the final coordinates in pdb format (boolean)
					save_format  : name of the extra binary file ( could be of the format: .dcd, .mdcrd ...) 
					rmsGradient  : root mean square gradient tolerance ( float )						
		'''
		_Optimizer = "ConjugatedGradient"
		_traj_name = None
		if "optmizer" 		 in _parameters: _Optimizer = _parameters["optmizer"]
		if "trajectory_name" in _parameters: _traj_name = _parameters["trajectory_name"]
		Gopt = GeometrySearcher(self.molecule,self.baseFolder,_trajName=_traj_name)		
		Gopt.ChangeDefaultParameters(_parameters)
		Gopt.Minimization(_Optimizer)
		Gopt.Finalize()
	#==================================================================
	def RelaxedSurfaceScan(self,_parameters,_plotParameters):
		'''
		Set up and execute one/two-dimensional relaxed surface scans 
		By the defualt the PKLs were saved on a child folder from the base path passed in the parameters, named "ScanTraj.ptGeo"
		The trajectory can be saved as files of the formats allowed by pDynamo 3.0
		Parameters:
			_parameters: python dict with parameters for simulation
				Mandatory :
					"ndim": 	number of reaction coordinates to be treated
					"ATOMS_RC1":list of atoms indices of the first reaction coordinate
					"nSteps_RC1": integer indicating the number of steps to scan for the first reaction coordinate
				Condirional:
					"ATOMS_RC2":list of atoms indices of the second reaction coordinate. Needed if "ndim = 2"
					"nSteps_RC2": integer indicating the number of steps to scan for the second reaction coordinate. Needed if "ndim = 2"
				Optinal   :
					"force_constant": Float indicating the constant value of energy penalty for the harmonic potential restriction function
					"force_constant_1" Specifies the force constant for the first reaction coordinate
					"force_constant_2" Specified the force constant for the second reaction coordinate
					"maxIterations": Number of maximum iteration for the geometry optimizations
					"optimizer": string containing the optimizer algorithm to be used in geometry optimization
					"dincre_RC1": float with the step increment for the first reaction coordinate ( Warning! If not passed, 0.0 will be assumed )
					"dincre_RC2": float with the step increment for the second reaction coordinate
					"MC_RC1": bool indicating whether to set mass constrained restrictions for the first reaction coordinate 
					"MC_RC2": bool indicating whether to set mass constrained restrictions for the second reaction coordinate 
					"rc_type_1": string containing the type for the first reaction coordinate ( Distance or Dihedral ) 
					"rc_type_2": string containing the type for the second reaction coordinate ( Distance or Dihedral )
					"adaptative": bool indicating wheter to activate or not the adaptative scheme for two-dimensional scans
					"save_format": format in which the trajectory will be saved, works only for 1D scans 
					"log_frequency": parameter for geometry optimization runs
			_plotParameters: python dict with paramters for post-analysis and plots
				All plot parameters are optionals. This dict can be passed as none, if so the plots will be perfomed with default parameters. 
					"contour_lines": integer indicating the number of contour lines to be used in two-dimensional plots
					"show":boolean indicating whether to display the plot before exiting. 
		'''
		#------------------------------------------------------------------
		#default varibales
		_Adaptative = False
		_Optmizer   = "ConjugatedGradient"
		MCR1 		= False
		MCR2 		= False
		rcType1     = "Distance"
		rcType2     = "Distance"
		nDims       = _parameters['ndim']
		dincre1     = 0.0
		dincre2     = 0.0
		nRC2        = 0
		#checking parameters
		if "dincre_RC1" in _parameters:	dincre1 	= _parameters["dincre_RC1"]
		if "dincre_RC2" in _parameters:	dincre2 	= _parameters["dincre_RC2"]	
		if "nSteps_RC2" in _parameters:	nRC2  		= _parameters["nSteps_RC2"]
		if "optmizer"   in _parameters:	_Optmizer   = _parameters["optmizer"]
		if "adaptative" in _parameters: _Adaptative = _parameters["adaptative"] 
		if "MC_RC1"     in _parameters: MCR1        = _parameters["MC_RC1"]
		if "MC_RC2"     in _parameters: MCR2        = _parameters["MC_RC2"]		
		if "rc_type_1"  in _parameters:	rcType1 	= _parameters["rc_type_1"]
		if "rc_type_2"  in _parameters:	rcType2 	= _parameters["rc_type_2"]
		
		#--------------------------------------------------------------------
		scan = SCAN(self.molecule,self.baseFolder,_Optmizer,ADAPTATIVE=_Adaptative)
		scan.ChangeDefaultParameters(_parameters)	
		#--------------------------------------------------------------------
		rc1 = ReactionCoordinate( _parameters["ATOMS_RC1"], MCR1, _type=rcType1 )
		rc1.SetInformation(self.molecule,dincre1)
		scan.SetReactionCoord(rc1)
		rc2 = None
		if nDims == 2:
			rc2 = ReactionCoordinate( _parameters["ATOMS_RC2"], MCR2, _type=rcType2 )
			rc2.SetInformation(self.molecule,dincre2)
			scan.SetReactionCoord(rc2)
			scan.Run2DScan(_parameters["nSteps_RC1"], _parameters["nSteps_RC2"] )
		else:
			scan.Run1DScan(_parameters["nSteps_RC1"])		
		scan.Finalize()		
		#================================================================
		#Set plor parameters
		cnt_lines 	= 12
		crd1_label	= rc1.label
		crd2_label	= ""
		show 		= False
		if nDims == 2: crd2_label = rc2.label
		#check parameters for plot
		if "contour_lines" 	in _plotParameters: cnt_lines = _plotParameters["contour_lines"]			
		if "show" 			in _plotParameters: show 	  = _plotParameters["show"] 
		#------------------------------------------------------------
		if 	 nDims 	== 2: TYPE = "2D"
		elif nDims 	== 1: TYPE = "1D"	
		#------------------------------------------------------------
		EA = EnergyAnalysis(_parameters['nSteps_RC1'],nRC2,_type=TYPE)
		EA.ReadLog( scan.baseName+"_scan{}D.log".format(nDims) ) 
		#-------------------------------------------------------------
		if 	 nDims == 2: EA.Plot2D(cnt_lines,crd1_label,crd2_label,show)
		elif nDims == 1: EA.Plot1D(crd1_label,show)		
	#=================================================================
	def MolecularDynamics(self,_parameters,_plotParameters):
		'''
		Set up and execute molecular dynamics simulations.
		Parameters:
			_parameters: python dict with parameters for simulation
				Mandatory:
					"MD_method": string containing the integrator algorithm name
					"protocol" : string indicating if is a normal run or for heating
					"nsteps"   : Number of steps to be taken in the simulation
				Optinal  :
					"temperature" 			  : float with the simulation temperature. If not passed we assume 300.15K as default.
					"coll_freq"  			  : integer with the colision frequency. Generally set for Langevin integrator. 
					"pressure"   			  : float with the simulation pressure. If not passed we assume 1.0bar as default.
					"pressure_coupling"		  : boolean indicating if is to control the simulation pressure.
					"temperature_scale_option": string with the type of temperature scaling. Default is 'linear' ( relevant for "heating" protocol)
					"temperature_scale"		  : float with the  temperature scaling step. Default is 10K  ( relevant for "heating" protocol)
					"start_temperatue"		  : float with the start temperature for heating protocol
					"timeStep"   			  : float indicating the size of integration time step. 0.001 ps is taken as default.					
					"sampling_factor"		  : integer indicating in which frequency to save/collect structure/data. default 0.
					"seed"					  : integer indicating the seed for rumdomness of the simulations.
					"log_frequency"     	  : integer indicating the frequency of the screen log output.
			_plotParameters: python dict with paramters for post-analysis and plots
				Optinal   :
					"show"					: whether to show the analysis plots in the simulation end.
					"calculate_distances"	: indicate if to calculate distances distributions of passed reaction coordinates
					"ATOMS_RC1"             : list of atoms for the first reaction coordinate to be analyzed 
					"ATOMS_RC2"             : list of atoms for the second reaction coordinate to be analyzed 
		'''		
		MDrun = MD(self.molecule,self.baseFolder,_parameters['MD_method'])		
		MDrun.ChangeDefaultParameters(_parameters)
		#---------------------------------------------------------------
		if "protocol" in _parameters:
			if   _parameters["protocol"] == "heating":  MDrun.HeatingSystem(_parameters['nsteps'],_parameters["sampling_factor"])
			elif _parameters["protocol"] == "sampling": MDrun.RunProduction(_parameters['nsteps'],_parameters["sampling_factor"])			
		#----------------------------------------------------------------
		if not _plotParameters == None:
			show = _plotParameters["show"]
			RCs  = None
			if "show" in _plotParameters: show = _plotParameters["show"]
			t_time = _parameters["nsteps"]*0.001
			DA = TrajectoryAnalysis(MDrun.trajectoryNameCurr,self.molecule,t_time)
			DA.CalculateRG_RMSD()
			DA.PlotRG_RMS(show)						
			if "calculate_distances" in _plotParameters:
				if _parameters["calculate_distances"] == True:
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
				Mandatory: 
					"ndim"                : integer indicating the number of restricted dimensions.
					"force_constant_1"    : float with the force constant applied for the harmonic potential in the first reaction coordinate
					"rc_type_1"           : string specifying the type of coordinate to be used for the first restricted dimension
					"rc_type_2"           : string specifying the type of coordinate to be used for the second restricted dimension
					"ATOMS_RC1"           : list with the atoms indices to set the first coordinate. 
					nsteps"               : number of moleculer dynamics steps to production run. 					
				Optinal :
					"ATOMS_RC2"           : list with the atoms indices to set the second coordinate. 
					"force_constant_2"    : float with the force constant applied for the harmonic potential in the second reaction coordinate
					"rc_type_2"           : string specifying the type of coordinate to be used for the second restricted dimension
					"temperature" 		  : float with the simulation temperature. If not passed we assume 300.15K as default.
					"coll_freq"  	      : integer with the colision frequency. Generally set for Langevin integrator. 
					"pressure"   		  : float with the simulation pressure. If not passed we assume 1.0bar as default.
					"pressure_coupling"	  : boolean indicating if is to control the simulation pressure.
					"timeStep"   		  : float indicating the size of integration time step. 0.001 ps is taken as default.
					"sampling_factor"     : integer with the save/sampling frequency of frames.
					"seed"				  : integer indicating the seed for rumdomness of the simulations.
					"log_frequency"       : integer indicating the frequency of the screen log output.
			_plotParameters:python dict with paramters for post-analysis and plots
				Optinal :
					"show"				  : whether to show the analysis plots in the simulation end.
					"calculate_distances" : whether to calculate and plot distribution analysis from the passed reaction coordinates										
		'''
		#----------------------------------------------------------------
		restraints = RestraintModel()
		self.molecule.DefineRestraintModel( restraints )		
		MCR1 = False
		MCR2 = False
		if "MC_RC1" in _parameters:	MCR1 = _parameters["MC_RC1"]
		if "MC_RC2" in _parameters:	MCR2 = _parameters["MC_RC2"]
		#--------------------
		rcType1 = "Distance"
		rcType2 = "Distance"
		if "rc_type_1" in _parameters: rcType1 = _parameters["rc_type_1"]
		if "rc_type_2" in _parameters: rcType2 = _parameters["rc_type_2"]
		#-------------------------------------------------------------------
		restrainDimensions = _parameters['ndim']
		forcK_1 = _parameters["force_constant_1"]
		#-------------------------------------------------------------------
		rc1 = ReactionCoordinate(_parameters["ATOMS_RC1"],MCR1,_type=rcType1)
		rc1.SetInformation(self.molecule,0)
		nDims = _parameters['ndim']
		rc2 = None
		if nDims == 2:
			rc2 = ReactionCoordinate(_parameters["ATOMS_RC2"],MCR2,_type=rcType2)
			rc2.SetInformation(self.molecule,0)
			forcK_2 = _parameters["force_constant_2"]
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
		MDrun.RunProduction(_parameters['nsteps'],_parameters["sampling_factor"],_Restricted=True)
		#-----------------------------------------------------------------		
		if not _plotParameters == None:
			t_time = _parameters["production_nsteps"]*MD.timeStep
			if "show" in _plotParameters: show = _plotParameters["show"]
			DA = TrajectoryAnalysis(MDrun.trajectoryNameCurr,self.molecule,t_time)
			DA.CalculateRG_RMSD()
			DA.PlotRG_RMS(show)				
			RCs = [rc1]
			if nDims > 1: RCs.append(rc2)							
			DA.DistancePlots(RCs,show)
			DA.ExtractFrames()
	#=======================================================================
	def UmbrellaSampling(self,_parameters,_plotParameters):
		'''
		Set up and execute umbrella sampling simulations and Free energy calculations for reaction path trajectory.
		Parameters:
			_parameters: python dict with parameters for simulation
				Mandatory:
				Optinal :
			_plotParameters:python dict with paramters for post-analysis and plots
				Mandatory:
				Optinal :
		'''
		#---------------------------------------
		MCR1 = False
		MCR2 = False		
		rcType1 = "Distance"
		rcType2 = "Distance"
		#---------------------------------------
		if "MC_RC1" in _parameters: 	MCR1 = True
		if "MC_RC2" in _parameters:		MCR2 = True
		#---------------------------------------
		_Restart 	= False
		_Adaptative = False
		if "restart" 	in _parameters:	_Restart 	= _parameters["restart"]
		if "adaptative" in _parameters:	_Adaptative = _parameters["adaptative"]
		#-------------------------------------------------------------------
		rc1 = ReactionCoordinate(_parameters["ATOMS_RC1"],MCR1,_type=rcType1)
		rc1.SetInformation(self.molecule,0)
		nDims = _parameters['ndim']
		rc2 = None
		if nDims == 2:
			rc2 = ReactionCoordinate(_parameters["ATOMS_RC2"],MCR2,_type=rcType2)
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
		#---------------
		USrun.Finalize()		
	#=========================================================================
	def PMFAnalysis(self,_parameters,_plotParameters):
		'''
		Calculate potential of mean force and Free energy from restricted molecular dynamics
		Parameters:
			_parameters: python dict with parameters for simulation
				Mandatory ( if not provided a "key-Error" will be thrown ):
				Optinal :
			_plotParameters:python dict with paramters for post-analysis and plots
				Mandatory ( if not provided a "key-Error" will be thrown ):
				Optinal :
		'''
		potmean = PMF( self.molecule, _parameters["source_folder"], self.baseFolder )
		potmean.CalculateWHAM(_parameters["xnbins"],_parameters["ynbins"],_parameters["temperature"])

		#================================================================
		#Set default plot parameters
		cnt_lines  = 12
		crd1_label = ""
		crd2_label = ""
		nRC2       = 0
		show       = False
		xwin       = 0
		ywin       = 0 
		#-----------------------------------------------------------------
		nDims = 1
		if _parameters["ynbins"] > 0:
			nDims = 2

		xlims = [ 0,  _parameters['xnbins'] ]
		ylims = [ 0,  _parameters['ynbins'] ]
		#-------------------------------------------------------------
		#check parameters for plot
		if "contour_lines" in _plotParameters:
			cnt_lines  = _plotParameters["contour_lines"]		
		if "xlim_list" in _plotParameters:
			xlims = _plotParameters["xlim_list"]
		if "ylim_list" in _plotParameters:
			ylims = _plotParameters["ylim_list"]
		if "ynbins" in _parameters:
			nRC2  = _parameters["ynbins"]
		if "show" in _plotParameters:
			show = True
		if "crd1_label" in _plotParameters:
			crd1_label = _plotParameters["crd1_label"]
		if "crd2_label" in _plotParameters:
			crd2_label = _plotParameters["crd2_label"]
		if "xwindows" in _plotParameters:
			xwin = _plotParameters["xwindows"]
		if "ywindows" in _plotParameters:
			ywin = _plotParameters["ywindows"]
		#------------------------------------------------------------
		if nDims == 2:
			TYPE = "WHAM2D"
		elif nDims == 1: 
			TYPE = "WHAM1D"	
		#------------------------------------------------------------
		# Plot PMF graphs
		EA = EnergyAnalysis(_parameters['xnbins'],nRC2,_type=TYPE)
		EA.ReadLog( potmean.baseName+".dat" ) 
		#-------------------------------------------------------------
		if nDims == 2:
			EA.Plot2D(cnt_lines,crd1_label,crd2_label,xlims,ylims,show)
		elif nDims == 1:
			EA.Plot1D(crd1_label,show)
		#-------------------------------------------
		#Plot Free energy of the calculated windows
		if nDims == 2:
			TYPE = "FE2D"
		elif nDims == 1: 
			TYPE = "FE1D"		
		#------------------------------------------
		EA = EnergyAnalysis(xwin,ywin,_type=TYPE)
		EA.ReadLog( potmean.baseName+"_FE.log" ) 
		#-------------------------------------------------------------
		if nDims == 2:
			EA.Plot2D(cnt_lines,crd1_label,crd2_label,xlims,ylims,show)
		elif nDims == 1:
			EA.Plot1D(crd1_label,show)
	#=========================================================================
	def NormalModes(self,_parameters):
		'''
		Simulation preset to calculate the normal modes and to write thr trajectory for a specific mode.
		Parameters:
			_parameters: python dict with parameters for simulation
				Mandatory ( if not provided a "key-Error" will be thrown ):
				Optinal :
			_plotParameters:python dict with paramters for post-analysis and plots
				Mandatory ( if not provided a "key-Error" will be thrown ):
				Optinal :
		'''	
		mode 		= 0
		temperature = 300.15
		Cycles 		= 10 
		Frames  	= 10 
		#-------------------------------
		if "temperature" in _parameters:
			temperature = _parameters["temperature"]
		if "cycles" in _parameters:
			Cycles = _parameters["cycles"]
		if "frames" in _parameters:
			Frames = _parameters["frames"]
		if "mode" in _parameters:
			mode   = _parameters["mode"]
		#-------------------------------
		NormalModes_SystemGeometry ( self.molecule, modify = ModifyOption.Project )
		if _mode > 0:
			trajectory = ExportTrajectory ( os.path.join (self.baseFolder, "NormalModes","ptGeo"), self.molecule )
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
				Mandatory ( if not provided a "key-Error" will be thrown ):
				Optinal :
			_plotParameters:python dict with paramters for post-analysis and plots
				Mandatory ( if not provided a "key-Error" will be thrown ):
				Optinal :
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
				Mandatory ( if not provided a "key-Error" will be thrown ):
				Optinal :
			_plotParameters:python dict with paramters for post-analysis and plots
				Mandatory ( if not provided a "key-Error" will be thrown ):
				Optinal :
		'''
		NEBrun = GeometrySearcher(self.molecule,self.baseFolder)

		#if there any parameters to be modified 
		NEBrun.ChangeDefaultParameters(_parameters)
		NEBrun.NudgedElasticBand(_parameters)


	#=========================================================================
	def SAW(self,_parameters):
		'''
		Set up and execute Self-Avoid Walking simulations to generate a reaction path trajectory 
		Parameters:
			_parameters: python dict with parameters for simulation
				Mandatory ( if not provided a "key-Error" will be thrown ):
				Optinal :
			_plotParameters:python dict with paramters for post-analysis and plots
				Mandatory ( if not provided a "key-Error" will be thrown ):
				Optinal :
		'''
	#=========================================================================	
	def TrajectoryPlots(self,_parameters,_plotParameters) :
		'''
		Produce graphical analysis from provided trajectories.
		Parameters:
			_parameters: python dict with parameters for simulation
				Mandatory ( if not provided a "key-Error" will be thrown ):
				Optinal :
			_plotParameters:python dict with paramters for post-analysis and plots
				Mandatory ( if not provided a "key-Error" will be thrown ):
				Optinal :
		'''

		pass
	#=========================================================================
	def EnergyPlots(self,_parameters,_plotParameters):
		'''
		Produce Energy plots from previus simulations log files
		Parameters:
			_parameters: python dict with parameters for simulation
				Mandatory ( if not provided a "key-Error" will be thrown ):
				Optinal :
			_plotParameters:python dict with paramters for post-analysis and plots
				Mandatory ( if not provided a "key-Error" will be thrown ):
				Optinal :
		'''		
		multiPlot = False
		ndim      = 1 
		crd1_label= "Reaction Coordinate #"
		crd2_label= "Reaction Coordinate #"
		cnt_lines = 0 
		ysize     = 0
		if "ysize" in _parameters: ysize = _parameters["ysize"]
		xlim      = [ 0, _parameters["xsize"] ]
		ylim 	  = [ 0, ysize ]
		show 	  = False
		#--------------------------------------------------------
		if "contour_lines" 	in _plotParameters: cnt_lines  = _plotParameters["contour_lines"]
		if "crd1_label" 	in _plotParameters: crd1_label = _plotParameters["crd1_label"]
		if "crd2_label" 	in _plotParameters:	crd2_label = _plotParameters["crd2_label"]
		if "xlim_list" 		in _plotParameters: xlim  	   = _plotParameters["xlim_list"]
		if "ylim_list" 		in _plotParameters: ylim       = _plotParameters["ylim_list"]
		if "show" 			in _plotParameters: show       = _plotParameters["show"]
		if "multiple_plot" 	in _plotParameters: multiPlot  = _plotParameters["multiple_plot"]		
		if ysize > 0: ndim = 2
		if "log_names" in _plotParameters: multiPlot = True 
		#--------------------------------------------------------
		EA = EnergyAnalysis(_parameters["xsize"],ysize,_type=_parameters["type"] )
		if multiPlot:
			for log in _parameters["log_names"]:
				EA.ReadLog( log )
				EA.MultPlot1D()
		else:
			EA.ReadLog( _parameters["log_name"] )
		#--------------------------------------------------------
		if 	 ndim == 1: EA.Plot1D(crd1_label,show)
		elif ndim == 2:	EA.Plot2D(cnt_lines,crd1_label,crd2_label,xlim,ylim,show)

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
