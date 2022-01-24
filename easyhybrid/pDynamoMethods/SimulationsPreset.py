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

#--------------------------------------------------------------
#loading pDynamo Libraries
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
		self.optmizer			= "ConjugatedGradient"
		self.samplingFactor 	= 1 # this is usually let to the default class value, unless the user want to modify
		self.nProcs 			= NmaxThreads
		#for restricted simulations
		
		#enviromental parameters and their default values 
		self.temperature = 300.15
		self.pressure 	 = 1

		#specific parameters for Energy refinement
		self.software 	= "pDynamo"
		self.methods  	= ["am1","rm1","pm3","pm6"] #some semiempirical methods as default for energy refinement

		
		#specific parameters for molecular dynamics run
		self.mdMethod 	= "Verlet"
		self.equiNsteps = 0
		self.prodNsteps = 0

		#specif parameters for the normal modes run
		self.NMcycles    = 10
		self.NMframes    = 20

		#--------------------------------------------------
		if not os.path.exists( self.baseFolder ):
			os.makedirs(self.baseFolder)

	#=======================================================================
	def Execute(self,_parameters):
		'''
		Function to call the class method to execute the preset simulation
		Also here is where some parameters default values can be updated 
		'''		
		#-------------------------------------------------------------
		if self.simulationType == "Energy_Refinement":			
			self.EnergyRefinement(_parameters)
		
		#-------------------------------------------------------------
		elif self.simulationType == "Geometry_Optimization":
			if "optmizer" in _parameters:
				self.optmizer = _parameters["optmizer"]

			self.GeometryOptimization(_parameters)

		#-------------------------------------------------------------
		elif self.simulationType == "Relaxed_Surface_Scan":			
			self.RelaxedSurfaceScan(_parameters)		

		#-------------------------------------------------------------
		elif self.simulationType == "Molecular_Dynamics":
			self.MolecularDynamics(_parameters)

		#-------------------------------------------------------------	
		elif self.simulationType == "Restricted_Molecular_Dynamics":			
			self.RestrictedMolecularDynamics(_parameters)

		#-------------------------------------------------------------
		elif self.simulationType == "Umbrella_Sampling":
			self.UmbrellaSampling(_parameters)
		
		#-------------------------------------------------------------
		elif self.simulationType == "Normal_Modes":
			if "temperature" in _parameters:
				self.temperature = _parameters['temperature']
			if "cycles" in _parameters:
				self.NMcycles = _parameters['cycles']
			if "frames" in _parameters:
				self.NMframes = _parameters['frames']	
			self.NormalModes()

		#-------------------------------------------------------------
		elif self.simulationType == "Delta_Free_Energy":
			if "temperature" in _parameters:
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

	

	#==================================================================
	def EnergyRefinement(self):
		'''
		Class method to set up and execute energy refinement using a series of methods
		'''

	#==================================================================
	def GeometryOptimization(self,_parameters):
		'''
		Class method to set up and execture the search of local minima for the system passed
		'''
		Gopt = GeometrySearcher(self.molecule,self.baseFolder)		
		Gopt.ChangeDefaultParameters(_parameters)
		Gopt.Minimization(self.optmizer)
		Gopt.Finalize()

	#==================================================================
	def RelaxedSurfaceScan(self,_parameters):
		'''
		Class method to set up and execute one/two-dimensional relaxed surface scans 
		'''
		#-------------------------------------------------------------------------
		scan = SCAN(self.molecule,self.baseFolder,self.optmizer,nprocs=self.nProcs)
		scan.ChangeDefaultParameters(_parameters)

		MCR1 = False
		MCR2 = False
		if "MC_RC1" in _parameters:
			MCR1 = True
		if "MC_RC2" in _parameters:
			MCR2 = True
		
		restraintDimensions = _parameters['ndim']

		scan.SetReactionCoord(_parameters['ATOMS_RC1'], _parameters['dincre_RC1'], MCR1)
		#--------------------------------------------------------------------------------
		if restraintDimensions == 2:
			scan.SetReactionCoord(_parameters['ATOMS_RC2'], _parameters['dincre_RC2'], MCR2)
			scan.RunTWODimensionSCAN(_parameters['nSteps_RC1'], _parameters['nSteps_RC2'] )
		else:
			scan.RunONEDimensionSCAN(_parameters['nSteps_RC1'])			
		#---------------------------------------------------------------------------------
		scan.Finalize()
		#.............

	#=================================================================
	def MolecularDynamics(self,_parameters):
		'''
		Class method to set up and execute molecular dynamics simulations.
		'''
		#-------------------------------------------------------------
		print(self.baseFolder)
		input()
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
		MDrun.Analysis()
		#...............
	#==================================================================
	def RestrictedMolecularDynamics(self,_parameters):
		'''
		Class method to set up and execute molecular dynamics simulations.
		#Ainda tenho que ver como functiona os arquivos de trajetórias na nova versão
		'''
		#----------------------------------------------------------------
		restraints = RestraintModel( )
		self.molecule.DefineRestraintModel( restraints )
		distances = []
		#----------------------------------------------------------------
		forcK = _parameters["forceC"]		
		restrainDimensions = _parameters['ndim']
		mdistance1 = False
		mdistance2 = False
		#---------------------------------------------------------------
		if len(_parameters["atoms"]) == 3:
			mdistance1 = True
		if len(_parameters["atoms"]) == 5:
			mdistance2 = True 
		if len(_parameters["atoms"]) == 6:
			mdistance1 = True
			mdistance2 = True 
		#--------------------------------------------------------------
		#task:
		#Improve the weights definitions to work with mass constraints
		weight1 = 1.0
		weight2 = -1.0 
		atoms = []
		
		for atom in range( len(_parameters["atoms"]) ):
			atoms.append( _parameters["atoms"][atom] )
		
		if restrainDimensions == 1:
			if mdistance1:
				distance = self.molecule.coordinates3.Distance(atoms[1],atoms[0]) - self.molecule.coordinates3.Distance(atoms[1],atoms[2]) 
				rmodel = RestraintEnergyModel.Harmonic( distance, forcK )
				restraint = RestraintMultipleDistance.WithOptions( energyModel = rmodel,  distances= [ [ atoms[1], atoms[0], weight1 ], [ atoms[1], atoms[2], weight2 ] ]) 
				restraints['ReactionCoord'] =  restraint
			else:
				distance = self.molecule.coordinates3.Distance(atoms[1],atoms[0])
				rmodel = RestraintEnergyModel.Harmonic( distance, forcK )
				restraint = RestraintDistance.WithOptions( energyModel = rmodel, point1= atoms[0], point2= atoms[1] )
				restraints['ReactionCoord'] =  restraint
		
		elif restrainDimensions == 2:
			if mdistance1 and mdistance2:
				distance1 = self.molecule.coordinates3.Distance( atoms[1],atoms[0] ) - self.molecule.coordinates3.Distance(atoms[1],atoms[2]) 
				rmodel1 = RestraintEnergyModel.Harmonic( distance1, forcK )
				restraint1 = RestraintMultipleDistance.WithOptions( energyModel = rmodel1,  distances= [ [ atoms[1], atoms[0], weight1 ], [ atoms[1], atoms[2], weight2 ] ]) 
				restraints['ReactionCoord1'] =  restraint1
				#
				distance2 = self.molecule.coordinates3.Distance(atoms[4],atoms[3]) - self.molecule.coordinates3.Distance(atoms[4],atoms[5]) 
				rmodel2 = RestraintEnergyModel.Harmonic( distance2, forcK )
				restraint2 = RestraintMultipleDistance.WithOptions( energyModel = rmodel2,  distances= [ [ atoms[4], atoms[3], weight1 ], [ atoms[4], atoms[5], weight2 ] ]) 
				restraints['ReactionCoord2'] =  restraint2
			elif mdistance1 and mdistance2 == False:

				distance1 = self.molecule.coordinates3.Distance( atoms[1],atoms[0] ) - self.molecule.coordinates3.Distance(atoms[1],atoms[2]) 
				rmodel1 = RestraintEnergyModel.Harmonic( distance1, forcK )
				restraint1 = RestraintMultipleDistance.WithOptions( energyModel = rmodel1,  distances= [ [ atoms[1], atoms[0], weight1 ], [ atoms[1], atoms[2], weight2 ] ]) 
				restraints['ReactionCoord1'] =  restraint1
				#
				distance2 = self.molecule.coordinates3.Distance(atoms[4],atoms[3])  
				rmodel2 = RestraintEnergyModel.Harmonic( distance2, forcK )
				restraint2 = RestraintDistance.WithOptions( energyModel = rmodel2, point1=atoms[3],point2=atoms[4] ) 
				restraints['ReactionCoord2'] =  restraint2
			else:
				distance1 = self.molecule.coordinates3.Distance( atoms[1],atoms[0] ) 
				rmodel1 = RestraintEnergyModel.Harmonic( distance1, forcK )
				restraint1 = RestraintDistance.WithOptions( energyModel = rmodel1, point1=atoms[0],point2=atoms[1] ) 
				restraints['ReactionCoord1'] =  restraint1
				#
				distance2 = self.molecule.coordinates3.Distance(atoms[4],atoms[3])  
				rmodel2 = RestraintEnergyModel.Harmonic( distance2, forcK )
				restraint2 = RestraintDistance.WithOptions( energyModel = rmodel2, point1=atoms[2],point2=atoms[3] ) 
				restraints['ReactionCoord2'] =  restraint2
				
		#----------------------------------------------------------------
		MDrun = MD(self.molecule,self.baseFolder,_parameters['MD_method'])
		MDrun.ChangeDefaultParameters(_parameters)
		MDrun.RunProductionRestricted(_parameters['equilibration_nsteps'],_parameters['production_nsteps'],200)
		MDrun.Analysis()
		#MDrun.DistAnalysis()

	#=======================================================================
	def UmbrellaSampling(self,_parameters):
		'''
		Class method to set up and execute umbrella sampling simulations and Free energy calculations for reaction path trajectory.
		#not finished method
		'''
		#---------------------------------------
		MCR1 = False
		MCR2 = False		
		if "MC_RC1" in _parameters:
			MCR1 = True
		if "MC_RC2" in _parameters:
			MCR2 = True
		#---------------------------------------
		USrun = US(self.molecule,self.baseFolder,_parameters['equilibration_nsteps'],_parameters['production_nsteps'],_parameters["MD_method"])
		USrun.ChangeDefaultParameters(_parameters)
		USrun.SetMode(_parameters["ATOMS_RC1"],MCR1)

		if _parameters["ndim"] == 1:
			USrun.RunONEDimensional(_parameters["trjFolder"],_parameters["samplingFactor"])
		elif _parameters["ndim"] == 2:
			USrun.SetMode(_parameters["ATOMS_RC2"],MCR2)
			USrun.RunTWODimensional(_parameters["trjFolder"],_parameters["samplingFactor"])
		
		prodfolders = self.baseFolder
		pmfRun = PMF(self.molecule,self.baseFolder)

	#=========================================================================
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

	#==========================================================================
	def DeltaFreeEnergy(self,_initCoord,_finalCoord):
		'''
		Calculate the free energy difference between two configurations of the system using the 
		statistical thermodynamics partition functions from through the normal modes calculations
		'''
		
		#initial Structure
		self.molecule.coordinates3 = ImportCoordinates3(_initCoord)
		e0 = self.molecule.Energy()
		NormalModes_SystemGeometry ( self.molecule, modify = ModifyOption.Project )
		Gibbs = [] 
		tdics = ThermodynamicsRRHO_SystemGeometry ( self.molecule 							,
                                                    pressure       = self.pressure       	,
                                                    symmetryNumber = self.symmetryNumber 	,
                                                    temperature    = self.temperature    	)
		Gibbs.append( tdics["Gibbs Free Energy"] )
    	# Final struct
		self.molecule.coordinates3 = ImportCoordinates3(_finalCoord)
		e1 = self.molecule.Energy()
		NormalModes_SystemGeometry ( self.molecule, modify = ModifyOption.Project )
    	 
		tdics = ThermodynamicsRRHO_SystemGeometry ( self.molecule 							,
                                                    pressure       = self.pressure       	,
                                                    symmetryNumber = self.symmetryNumber 	,
                                                    temperature    = self.temperature    	)
		Gibbs.append( tdics["Gibbs Free Energy"] )

    	#Analyze the results 


	#=========================================================================	
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

	#=========================================================================
	def SAW(self):
		'''
		Class method to set up and execute Self-Avoid Walking simulations to generate a reaction path trajectory 
		'''
	#========================================================================
	def SimulatingAnnealing(self):
		'''
		'''
		pass
	#========================================================================
	def SMD(self):
		'''
		'''
		pass
	#;;;;;;;;;
#=============================================================================
#========================END OF THE FILE======================================
#=============================================================================