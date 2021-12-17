#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = MolecularDynamics.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#==============================================================================

#---------------------------------------
#importing libraries
import os
import time
#import sys

#importing our library functions
import commonFunctions

# pDynamo
import pBabe
import pCore
import pMolecule
import pMoleculeScripts
#---------------------------------------

#==============================================================================
class MD:
    '''
    Class to set up Molecular Dynamics Sumulations.
    '''
    #.---------------------------------------
    def __init__(self,_system,_baseFolder,_method):
        '''
        Default constructor. 
        Receives a list of paramters to set the simulations.
        '''        
        #Important parameters that are recurrently wanted to be change by the user
        self.molecule               = _system
        self.trajectoryNameEqui     = os.path.join(_baseFolder,"equilibration","trj")
        self.trajectoryNameProd     = os.path.join(_baseFolder,"production","trj")
        self.trajectoryNameSoft     = os.path.join(_baseFolder,"soft_constr","trj")
        self.trajectoryNameCurr     = self.trajectoryNameEqui
        self.prodNsteps             = 20000 # seting the default for umbrella sampling
        self.equiNsteps             = 5000  # seting the default for umbrella sampling
        self.timeStep               = 0.001
        self.temperature            = 300.15
        self.pressureControl        = False
        self.algorithm              = _method
        self.samplingFactor         = 20
        self.logFreq                = 200
        self.seed                   = 304434242
        self.logname                = ""
        self.softConstraint         = False # boolean flags signlizing whether the system has soft constraints or not              
        #Parameters to be used internally for information
        self.state                  = "non-setted"  # string sinalizando o status da instancia dessa classe. 
        self.equilibrated           = False         # Signlizes the whether the system passed through a succesfull equilibrtion simulation 
        self.isQuantum              = False         # Boolean flag signilizing if some extent of the system is treated with a quatum energy model
        self.equi_time              = self.equiNsteps * self.timeStep            # default time for equilibration in ps
        self.prod_time              = self.prodnsteps * self.timeStep # production time where the data must to be collected
        self.outputDCD              = True          # Boolean flag signalizing whether the saved trajectory must to be save as a ".dcd" file.
        #Default constants less acessible by the users
        self.collFreq               = 25.0 
        self.pressureCoupling       = 2000.0       
        self.temperatureScaleFreq   = 0.1
        self.pressure               = 1.0   
        self.temperatureScaleOption = "linear"
        self.startTemperature       = 10      
        #Parameters that can be usefull to change for umbrella sampling
        self.maxItQC                = 250
        self.energyTolQC            = 2.0e-4
        self.densityTol             = 2.0e-8    
        #Setting parameters based on information that we collected on the instance construction
        self.RNG                    = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( self.seed ) )
    
    #.---------------------------------------    
    def ChangeDefaultParameters(self,_parameters):
        '''
        Class method to set more specifc parameters.
        
        For Bachega to think:
            there are some of those parameters that the commom user will not have interest in change
            and then you may want to not acces such method through EasyHybrid GUI.
            
        '''        
        if 'coll_freq'                  in _parameters:
            self.collFreq               = _parameters['coll_freq']
        if 'pressure_coupling'          in _parameters:
            self.pressureCoupling       = _parameters['pressure_coupling'] 
        if 'temperature_scale'          in _parameters:
            self.temperatureScaleFreq   = _parameters['temperature_scale']
        if 'maxIterations_QC'           in _parameters:
            self.maxItQC                = _parameters['maxIterations_QC']
        if 'energy_tolerance'           in _parameters:
            self.energyTolQC            = _parameters['energy_tolerance']
        if 'density_tolerancen'         in _parameters:
            self.densityTol             = _parameters['density_tolerance']
        if 'timeStep'                   in _parameters:
            self.timeStep               = _parameters['timeStep']
        if 'pressureControl'            in _parameters: # the presence of the ket activates this boolean flag
            self.pressureControl        = True
        if 'save_frequency'             in _parameters: 
            self.samplingFactor         = _parameters['saveFreq']
        if 'log_frequency'              in _parameters:
            self.logFrequency           = _parameters['log_frequency']
        if 'temperature_scale_option'   in _parameters:
            self.temperatureScaleOption = _parameters['temperature_scale_option']

    #.---------------------------------------    
    def HeatingSystem(self,_nsteps):
        '''
        Run a Velocity Verlet molecular dynamics simulation to gradually 
        make the system reach certain temperature. 
        '''

        Log = TextLogFileWriter(self.logname)

        VelocityVerletDynamics_SystemGeometry(  self.molecule                                        ,
                                                log = Log                                            ,
                                                logFrequency              = self.logFreq             ,
                                                normalDeviateGenerator    = self.RNG                 ,
                                                steps                     = _nsteps                  ,
                                                timeStep                  = self.timeStep            ,
                                                temperatureScaleFrequency = self.temperatureScaleFreq,
                                                temperatureScaleOption    = self.temperatureScale    ,
                                                temperatureStart          = self.startTemperature    ,
                                                temperatureStop           = self.temperature         )
    
        self.state = "Heated"
    #.---------------------------------------
    def RunEquilibration(self,_equiSteps):
        '''
        Run a molecular dynamics simulation for equilibration of the system
        '''

        self.nsteps             = _equiSteps
        self.trajectoryNameCurr = self.trajectoryNameEqui
        self.logname            = os.path.join( self.trajectoryNameCurr, "equi.log" )

        if self.algorithm == "Verlet":
            self.runVerlet()
        elif self.algorithm == "LeapFrog":
            self.runLeapFrog()
        elif self.algorithm == "Langevin":
            self.runLangevin()

        if self.CheckSimulation():
            self.state = "Equilibrated"
   
    
    #.---------------------------------------    
    def RunProduction(self,_prodSteps,):
        '''
        Run a molecular dynamics simulation for data collection.
        '''
        self.nsteps             = _prodSteps
        self.trajectoryNameCurr = self.trajectoryNameProd         
        self.logname            = os.path.join( self.trajectoryNameCurr, "prod.log" )

        
        if self.algorithm == "Verlet":
            self.runVerlet()
        elif self.algorithm == "LeapFrog":
            self.runLeapFrog()
        elif self.algorithm == "Langevin":
            self.runLangevin()

        if self.CheckSimulation("production"):
            self.state = "Sampled"

        if self.outputDCD:
            DCDTrajectory_FromSystemGeometryTrajectory(self.trajectoryNameCurr+".dcd", self.trajectoryNameCurr , self.molecule )
  
   
    #.---------------------------------------
    def RunProductionRestricted(self,_equiSteps,_prodSteps,_samplingFactor):
        '''
        Run a simulation with the system having soft constrains defined.
        Designed for Umbrella Sampling and Steered Molecular Dynamics routines.
        '''
        self.nsteps             = _prodSteps
        self.samplingFactor     = _samplingFactor
        self.trajectoryNameCurr = self.trajectoryNameProd         
        self.logname            = os.path.join( self.trajectoryNameCurr, "prod.log" )

        if self.outputDCD:
            DCDTrajectory_FromSystemGeometryTrajectory(  self.trajectoryNameCurr+".dcd", self.trajectoryNameCurr , self.molecule )
 
    
    #.---------------------------------------
    def CheckSimulation(self,_type):
        '''
        Class method to check whether the molecular dynamics ended well. 
        In the case of the equilibration runs, check whether the simulations reached satisfactorly conditions. 

        Tasks: 
            1. Develop algorithm that decides if the simulation go well
            2. attach plotting functionalities for visualization 
            3. write the method
        '''
        pass   

    #.---------------------------------------
    def runVerlet(self):
        '''
        Execute velocity verlet molecular dynamics from pDynamo methods. 
        '''
        trajectory          = ExportTrajectory( self.trajectoryNameCurr, self.molecule ) 
        trajectorySC        = None       
        Log                 = TextLogFileWriter(self.logname)
        trajectory_lists    = [] 

        if self.softConstraint:
            trajectorySC = ExportTrajectory(self.trajectoryNameSoft, self.molecule)

        VelocityVerletDynamics_SystemGeometry(self.system                                                   ,
                                log                         = Log                                           , 
                                logFrequency                = self.logFreq                                  ,
                                normalDeviateGenerator      = self.RNG                                      ,
                                steps                       = self.nsteps                                   ,
                                timeStep                    = self.timeStep                                 ,
                                temperatureScaleFrequency   = self.temperatureScaleFreq                     ,
                                temperatureScaleOption      = "constant"                                    ,
                                trajectories                = [ ( trajectory, self.samplingFactor ) ]       ,
                                temperatureStart            =   self.temperature                            )

        #.---------------------------------------
    def runLeapFrog(self):
        '''
        Execute Leap Frog molecular dynamics from pDynamo methods.
        '''

        trajectory  = ExportTrajectory(self.trajectoryNameCurr, self.molecule)
        Log         = TextLogFileWriter(self.logname)

        LeapFrogDynamics_SystemGeometry(self.system                                                         ,
                                        trajectories            = [ ( trajectory, self.samplingFactor ) ]   ,
                                        log                     = log                                       ,
                                        logFrequency            = self.logFreq                              ,
                                        normalDeviateGenerator  = self.RNG                                  ,
                                        pressure                = self.pressure                             ,
                                        pressureCoupling        = self.pressureCoupling                     ,
                                        steps                   = self.nsteps                               ,
                                        timeStep                = self.timeStep                             ,
                                        temperature             = self.temperature                          , 
                                        temperatureCoupling     = 0.1                                       ) 
                                
        
        if self.outputDCD:
            DCDTrajectory_FromSystemGeometryTrajectory(  self.trajectoryNameCurr+".dcd", self.trajectoryNameCurr , self.molecule )
        
    
    #.---------------------------------------    
    def runLangevin(self):
        '''
        Execute Langevin molecular dynamics from pDynamo methods.
        '''
        trajectory  = ExportTrajectory(self.trajectoryNameCurr, self.molecule)
        Log         = TextLogFileWriter(self.logname)
        
        LangevinDynamics_SystemGeometry ( self.system                           ,
                                          collisionFrequency     = self.collFreq,
                                          logFrequency           = self.logFreq ,
                                          normalDeviateGenerator = self. ,
                                          steps                  = _NSteps ,
                                          temperature            =   300.0 ,
                                          timeStep               =   0.001 ,
                                          trajectories = [ ( trajectory, _NSave ) ] )


        if self.outputDCD:
            DCDTrajectory_FromSystemGeometryTrajectory(self.trajectoryNameCurr+".dcd",self.trajectoryNameCurr,self.molecule) 

    #.-----------------------------------------
    def UnitTest(self):
        pass

#. End of class MD
#==============================================
       
    
    