#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = SimulationsPreset.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

import os, glob
import pCore
import pMolecule
import pBabel


#=============================================================
class Simulation:
	'''
	Class to set up preset simulations to be perfomed
	'''
	def __init__(self,_system,_simulationType):
		'''
		'''
		self.molecule 			= _system
		self.simulationType 	= _simulationType
		self.wall_time 			= 0 
		self.baseFolder  		=  os.getcwd() # the baseFolder for the simulations will be the current dir for now.
		self.MAXnprocs 			= 1 # maximum number of virtual threads to be used in the simulations
		self.coorddinatesFolder = "" # Name of the folder containing the pkls to be read. Used in more than one preset here
		self.logFreq 			= 1
		self.optmizer			= "ConjugatedGradient"
		self.samplingFactor 	= 1 # this is usually let to the default class value, unless the user want to modify
		#for restricted simulations
		self.restrainDimensions = 1
		
		#enviromental parameters and their default values 
		self.temperature = 300.15
		self.pressure 	 = 1

		#specific parameters for Energy refinement
		self.software 	= "pDynamo"
		self.methods  	= ["am1","rm1","pm3","pm6"] #some semiempirical methods as default for energy refinement

		
		#specific parameters for molecular dynamics run
		self.mdMethod 	= "Verlet"
		self.equiNsteps = 5000
		self.prodNsteps = 20000 

		#specif parameters for the normal modes run
		self.NMcycles    = 10
		self.NMframes    = 20


	#-------------------------------------------------------
	def Execute(self,_parameters):
		'''
		Function to call the class method to execute the preset simulation
		Also here is where some parameters default values can be updated 
		'''
		
		#-------------------------------------------------------------
		elif self.simulationType == "Energy_Refinement":
			if "methods" in _parameters:
				self.methods = _parameters
			if "software" in _parameters:
				self.software = _ṕarameters:
			if "nprocs" in  _parameters:
				self.MAXnprocs = _parameters['nprocs']
			if "coorddinates_folder" in _parameters:
				self.coorddinatesFolder = _parameters["coorddinates_folder"]
			self.EnergyRefinement()
		
		#-------------------------------------------------------------
		elif self.simulationType == "Geometry_Optimization":
			if "optmizer" in _parameters:
				self.optmizer = _parameters["optmizer"]

			if len(_parameters) > 1:
				self.GeometryOptimization(_parameters)
			else 
				self.GeometryOptimization(None)

		#-------------------------------------------------------------
		elif self.simulationType == "Relaxed_Surface_Scan":			
			self.RelaxedSurfaceScan(_parameters)		

		#-------------------------------------------------------------
		elif self.simulationType == "Molecular_Dynamics":
			if "MD_method" in _parameters:
				self.mdMethod = _parameters['MD_method'] 
			if "temperature" in _parameters:
				self.temperature = _parameters['temperature']
			if "production_nsteps" in _parameters:
				self.prodSteps 	=_parameters['production_nsteps']
			if "equilibration_nsteps" in _parameters:
				self.equiNsteps = _parameters['equilibration_nsteps']

			self.MolecularDynamics(_parameters)

		#-------------------------------------------------------------	
		elif self.simulationType == "Restricted_Molecular_Dynamics":
			self.RestrictedMolecularDynamics()

		#-------------------------------------------------------------
		elif self.simulationType == "Umbrella_Sampling":
			if "production_nsteps" in _parameters:
				self.prodNsteps = _parameters['production_nsteps']
			if "equilibration_nsteps" in _parameters:
				self.equiNsteps = _parameters['equilibration_nsteps']
			if "temperature" in _parameters
				self.temperature = _parameters['temperature']
			if "MD_method" in _parameters:
				self.mdMethod = _parameters['MD_method'] 
			
			self.UmbrellaSampling(_parameters)
		
		#-------------------------------------------------------------
		elif self.simulationType == "Normal_Modes":
			if "temperature" in _parameters
				self.temperature = _parameters['temperature']
			if "cycles" in _parameters:
				self.NMcycles = _parameters['cycles']
			if "frames" in _parameters:
				self.NMframes = _parameters['frames']	
			self.NormalModes()

		#-------------------------------------------------------------
		elif self.simulationType == "Delta_Free_Energy":
			if "temperature" in _parameters
				self.temperature = _parameters['temperature']
			_ic = _parameters['initial_coordinates']
			_fc = _parameters['final_coordinates']
			self.DeltaFreeEnergy(_ic,_fc)

		#-------------------------------------------------------------
		elif self.simulationType == "NEB":
			self.NEB()

		#-------------------------------------------------------------
		elif self.simulationType == "SAW":
			self.SAW()

		#-------------------------------------------------------------
		elif self.simulationType == "Simulating_Annealing":
			self.SimulatingAnnealing()

		#-------------------------------------------------------------
		elif self.simulationType == "Steered_Molecular_Dynamics":
			self.SMD()

	

	#-------------------------------------------------------
	def EnergyRefinement(self):
		'''
		Class method to set up and execute energy refinement using a series of methods
		'''

	#-------------------------------------------------------
	def GeometryOptimization(self,_parameters):
		'''
		Class method to set up and execture the search of local minima for the system passed
		'''
		Gopt = GeometrySearcher(self.molecule,self.baseFolder)
		
		#If there more parameters passed for this type of simulation the code understands that the user wants to modufy some deafult values
		if _parameters not None:
			Gopt.ChengeDefaultParameters(_parameters)

		Gopt.Minimization(self.optmizer)


	#_------------------------------------------------------
	def RelaxedSurfaceScan(self,_parameters):
		'''
		Class method to set up and execute one/two-dimensional relaxed surface scans 
		'''
		scan = SCAN(self.molecule,self.baseFolder,self.optmizer)

		if _parameters['change_parameters'] == "change":
			scan.ChengeDefaultParameters(_parameters)

		MCR1 = False
		MCR2 = False
		if _parameters["Mass_Constraint_RC1"]
			MCR1 = True
		if _parameters["Mass_Constraint_RC2"]
			MCR2 = True

		scan.SetReactionCoord(_parameters['ATOMS_RC1'], _parameters['Distance_Step_RC1'], MCR1)
		
		if self.restrainDimensions == 2:
			scan.SetReactionCoord(_parameters['ATOMS_RC2'], _parameters['Distance_Step_RC2'], MCR2)
			scan.RunTwoDimensionalSCAN(_parameters['nSteps_RC1'], _parameters['nSteps_RC2'] )
		else:
			scan.RunTwoDimensionalSCAN(_parameters['nSteps_RC1'])


	#_------------------------------------------------------
	def MolecularDynamics(self,_parameters):
		'''
		Class method to set up and execute molecular dynamics simulations.
		'''
		MDrun = MD(self.molecule,self.baseFolder,self.mdMethod)
		
		MDrun.temperature 	= self.temperature
		MDrun.pressure 		= self.pressure

		#If there is some key in _parameters set to modify some varibles then this is done here
		#If there is any this next command will have no effect
		MDrun.ChengeDefaultParameters(_parameters)

		if "protocol" in _parameters:
			if _parameters["protocol"] 		== "heating":
				MDrun.HeatingSystem(self.prodnsteps)
			elif _parameters["protocol"] 	== "equilibration":
				MDrun.RunEquilibration(self.equiNsteps)
			elif _parameters["protocol"] 	== "production":
				MDrun.RunEquilibration(self.equiNsteps)
				MDrun.RunProduction(self.prodnsteps)


	#_------------------------------------------------------
	def RestrictedMolecularDynamics(self):
		'''
		Class method to set up and execute molecular dynamics simulations.
		#Ainda tenho que ver como functiona os arquivos de trajetórias na nova versão
		'''

	#_------------------------------------------------------
	def UmbrellaSampling(self):
		'''
		Class method to set up and execute umbrella sampling simulations and Free energy calculations for reaction path trajectory.
		'''
		USrun = US(self.molecule,self.molecule,self.equiNsteps,self.prodNsteps,self.mdMethod)

		MCR1 = False
		MCR2 = False
		
		if _parameters["Mass_Constraint_RC1"]
			MCR1 = True
		if _parameters["Mass_Constraint_RC2"]
			MCR2 = True

		if self.restrainDimensions == 2:

	#_------------------------------------------------------
	def NormalModes(self,_mode):
		'''
		Simulation preset to calculate the normal modes and to write thr trajectory for a specific mode.
		'''	
		NormalModes_SystemGeometry ( self.molecule, modify = ModifyOption.Project )
		if _mode > 0:
			trajectory = ExportTrajectory ( os.path.join (self.baseFolder, "NormalModes","trj"), self.molecule )
			NormalModesTrajectory_SystemGeometry(	self.molecule					,
                                       			 	trajectory          			,
                                       				mode        = _mode				,
                                       				cycles      = self.NMcycles    	,
                                       				frames      = self.NMframes 	,
                                       				temperature = self.temperature )

	#_------------------------------------------------------
	def DeltaFreeEnergy(self,_initCoord,_finalCoord):
		'''
		Calculate the free energy difference between two configurations of the system using the 
		statistical thermodynamics partition functions from through the normal modes calculations
		'''
		
		#initial Structure
		self.molecule = Unpickle(_initCoord):
		e0 = self.molecule.Energy()
    	NormalModes_SystemGeometry ( self.molecule, modify = ModifyOption.Project )
    	Gibbs = [] 
		tdics = ThermodynamicsRRHO_SystemGeometry ( self.molecule 							,
                                                    pressure       = self.pressure       	,
                                                    symmetryNumber = self.symmetryNumber 	,
                                                    temperature    = self.temperature    	)
    	Gibbs.append( tdics["Gibbs Free Energy"] )
    	# Final struct
    	self.molecule = Unpickle(_finalCoord):
		e1 = self.molecule.Energy()
    	NormalModes_SystemGeometry ( self.molecule, modify = ModifyOption.Project )
    	 
		tdics = ThermodynamicsRRHO_SystemGeometry ( self.molecule 							,
                                                    pressure       = self.pressure       	,
                                                    symmetryNumber = self.symmetryNumber 	,
                                                    temperature    = self.temperature    	)
    	Gibbs.append( tdics["Gibbs Free Energy"] )

    	#Analyze the results 


	#_------------------------------------------------------	
	def NEB(self,_parameters):
		'''
		Class method to set up and execute Nudget Elastic Band simulations to generate a reaction path trajectory 
		'''
		NEBrun = GeometrySearcher(self.system,self.baseFolder)

		#if there any parameters to be modified 
		NEBrun.ChengeDefaultParameters(_parameters)

		NEBrun.NudgedElasticBand( 	_parameters['initial_coordinates'], 	\
									_parameters['final_coordinates'] ,  	\
									_parameters['NEB_nbins'],				\
									_parameters["RMS_growing_intial_string"] 
								)

	#_------------------------------------------------------
	def SAW(self):
		'''
		Class method to set up and execute Self-Avoid Walking simulations to generate a reaction path trajectory 
		'''
	#------------------------------------------------------
	def SimulatingAnnealing(self):
		'''
		'''
		pass
	#_------------------------------------------------------
	def SMD(self):
		'''
		'''
		pass
