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
#import sys
#----------------------------------------
# pDynamo
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
#---------------------------------------

#**************************************************************************
class MD:
    '''
    Class to set up Molecular Dynamics Sumulations.
    '''
    #.---------------------------------------
    def __init__(self,_system,_baseFolder,_integrator):
        '''
        Default constructor. 
        Receives a list of paramters to set the simulations.
        '''        
        #Important parameters that are recurrently wanted to be change by the user
        self.molecule               = _system
        self.baseName               = _baseFolder       
        self.trajectoryNameProd     = os.path.join(_baseFolder,"trajectory.ptGeo")
        self.trajectoryNameSoft     = _baseFolder+"trajectory.ptRes" # this naming scheme makes sense for the umbrella sampling runs
        self.trajectoryNameCurr     = self.trajectoryNameProd
        self.algorithm              = _integrator
        self.saveFormat             = None # binary file format to save the trajectory
        self.Nsteps                 = 20000 
        self.timeStep               = 0.001
        self.temperature            = 300.15
        self.pressureControl        = False
        self.samplingFactor         = 0
        self.logFreq                = 200
        self.seed                   = 3029202042
        self.softConstraint         = False # boolean flags signlizing whether the system has soft constraints or not              
        #Default constants less acessible by the users
        self.collFreq               = 25.0 
        self.pressureCoupling       = 2000.0       
        self.temperatureScaleFreq   = 10
        self.pressure               = 1.0   
        self.temperatureScaleOption = "linear"
        self.startTemperature       = 10               
        #Setting parameters based on information that we collected on the instance construction
        self.RNG                    = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( self.seed ) )
        if not os.path.exists(_baseFolder): os.makedirs(_baseFolder)
    #===============================================================================    
    def ChangeDefaultParameters(self,_parameters):
        '''
        Class method to set more specifc parameters.           
        '''     
        if "temperature"                in _parameters: self.temperature            = _parameters["temperature"]   
        if "start_temperature"          in _parameters: self.startTemperature       = _parameters["start_temperature"]
        if "coll_freq"                  in _parameters: self.collFreq               = _parameters["coll_freq"]
        if "pressure"                   in _parameters: self.pressure               = _parameters["pressure"]
        if "pressure_coupling"          in _parameters: self.pressureCoupling       = _parameters["pressure_coupling"] 
        if "temperature_scale"          in _parameters: self.temperatureScaleFreq   = _parameters["temperature_scale"]
        if "timeStep"                   in _parameters: self.timeStep               = _parameters["timeStep"]
        if "sampling_factor"            in _parameters: self.samplingFactor         = _parameters["sampling_factor"]
        if "log_frequency"              in _parameters: self.logFrequency           = _parameters["log_frequency"]
        if "temperature_scale_option"   in _parameters: self.temperatureScaleOption = _parameters["temperature_scale_option"]
        if "seed"                       in _parameters:
            self.seed = _parameters["seed"]
            self.RNG  = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( self.seed ) )
    #=============================================================================================    
    def HeatingSystem(self,_nsteps):
        '''
        Run a Velocity Verlet molecular dynamics simulation to gradually 
        make the system reach certain temperature. 
        '''
        self.nsteps             = _nsteps
        self.trajectoryNameProd = os.path.join(self.baseName,"heating.ptGeo")  
        self.trajectoryNameCurr = self.trajectoryNameProd 
        trajectory              = ExportTrajectory( self.trajectoryNameCurr, self.molecule,log=None )

        if not os.path.exists( self.trajectoryNameCurr ): os.makedirs( self.trajectoryNameCurr )
        #---------------------------------------------------------------------------------------------
        VelocityVerletDynamics_SystemGeometry(self.molecule                                                 ,
                                              logFrequency              = self.logFreq                        ,
                                              normalDeviateGenerator    = self.RNG                            ,
                                              steps                     = self.nsteps                         ,
                                              timeStep                  = self.timeStep                       ,
                                              trajectories              = [(trajectory,self.samplingFactor)]  ,
                                              temperatureScaleFrequency = self.temperatureScaleFreq           ,
                                              temperatureScaleOption    = self.temperatureScaleOption         ,
                                              temperatureStart          = self.startTemperature               ,
                                              temperatureStop           = self.temperature                    )            
    #===============================================================================================    
    def RunProduction(self,_prodSteps,_samplingFactor,_Restricted=False):
        '''
        Run a molecular dynamics simulation for data collection.
        '''
        self.softConstraint     = _Restricted
        self.nsteps             = _prodSteps
        self.trajectoryNameCurr = self.trajectoryNameProd         
        self.samplingFactor     = _samplingFactor
        if not os.path.exists( self.trajectoryNameCurr ): os.makedirs( self.trajectoryNameCurr )
        if   self.algorithm == "Verlet":      self.runVerlet()
        elif self.algorithm == "LeapFrog":    self.runLeapFrog()
        elif self.algorithm == "Langevin":    self.runLangevin()       
    #===================================================================================================
    def runVerlet(self):
        '''
        Execute velocity verlet molecular dynamics from pDynamo methods. 
        '''
        trajectory      = ExportTrajectory( self.trajectoryNameCurr, self.molecule,log=None )         
        trajectory_list = []

        if self.softConstraint:
            trajSoft = ExportTrajectory(self.trajectoryNameSoft, self.molecule,log=None)
            trajectory_list = [ ( trajectory, self.samplingFactor ), (trajSoft, 1) ]
        else:
            trajectory_list = [ ( trajectory, self.samplingFactor ) ]
        
        VelocityVerletDynamics_SystemGeometry(  self.molecule                                           ,
                                                logFrequency                = self.logFreq              ,
                                                normalDeviateGenerator      = self.RNG                  ,
                                                steps                       = self.nsteps               ,
                                                timeStep                    = self.timeStep             ,
                                                temperatureScaleFrequency   = self.temperatureScaleFreq ,
                                                temperatureScaleOption      = "constant"                ,
                                                trajectories                = trajectory_list           ,
                                                temperatureStart            = self.temperature          )

    #====================================================================================================
    def runLeapFrog(self):
        '''
        Execute Leap Frog molecular dynamics from pDynamo methods.
        '''
        #--------------------------------------------------------------------------------
        trajectory  = ExportTrajectory(self.trajectoryNameCurr, self.molecule,log=None)       
        trajectory_list = []
        #--------------------------------------------------------------------------------
        if self.softConstraint:
            trajSoft = ExportTrajectory(self.trajectoryNameSoft, self.molecule)
            trajectory_list = [ ( trajectory, self.samplingFactor ), ( trajSoft, 1 ) ]
        else:
            trajectory_list = [ ( trajectory, self.samplingFactor ) ]
        #--------------------------------------------------------------------------------
        LeapFrogDynamics_SystemGeometry(self.molecule                                   ,
                                        trajectories            = trajectory_list       ,
                                        logFrequency            = self.logFreq          ,
                                        normalDeviateGenerator  = self.RNG              ,
                                        pressure                = self.pressure         ,
                                        pressureCoupling        = self.pressureCoupling ,
                                        steps                   = self.nsteps           ,
                                        timeStep                = self.timeStep         ,
                                        temperatureControl      = True                  ,
                                        temperature             = self.temperature      , 
                                        temperatureCoupling     = 0.1                   ) 

    
    #======================================================================================
    def runLangevin(self):
        '''
        Execute Langevin molecular dynamics from pDynamo methods.
        '''
        #-----------------------------------------------------------------------------
        trajectory  = ExportTrajectory(self.trajectoryNameCurr, self.molecule)
        trajectory_list = []
        #-----------------------------------------------------------------------------
        if self.softConstraint:
            trajSoft = ExportTrajectory(self.trajectoryNameSoft, self.molecule)
            trajectory_list = [ ( trajectory, self.samplingFactor ), ( trajSoft, 1 ) ]
        else:
            trajectory_list = [ ( trajectory, self.samplingFactor ) ]
        #-----------------------------------------------------------------------------
        LangevinDynamics_SystemGeometry ( self.molecule                             ,
                                          collisionFrequency     = self.collFreq    ,
                                          logFrequency           = self.logFreq     ,
                                          normalDeviateGenerator = self.RNG         ,
                                          steps                  = self.nsteps      ,
                                          temperature            = self.temperature ,
                                          timeStep               = self.timeStep    ,
                                          trajectories           = trajectory_list  )

    #=====================================================================================
         

#===============================================================================#
#. End of class MD =============================================================#
#===============================================================================#
       
    
    
