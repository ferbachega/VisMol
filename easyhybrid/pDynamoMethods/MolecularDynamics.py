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
#importing our library functions
import commonFunctions
from LogFile import LogFile
from ReactionCoordinate import *
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
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import seaborn as sns
#---------------------------------------

#**************************************************************************
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
        self.baseName               = _baseFolder       
        self.trajectoryNameEqui     = os.path.join(_baseFolder,"equilibration.ptGeo")
        self.trajectoryNameProd     = os.path.join(_baseFolder,"production.ptGeo")
        self.trajectoryNameSoft     = _baseFolder+"restricted.ptRes"
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
        self.softConstraint         = False # boolean flags signlizing whether the system has soft constraints or not              
        self.outputDCD              = True          # Boolean flag signalizing whether the saved trajectory must to be save as a ".dcd" file.
        #Default constants less acessible by the users
        self.collFreq               = 25.0 
        self.pressureCoupling       = 2000.0       
        self.temperatureScaleFreq   = 10
        self.pressure               = 1.0   
        self.temperatureScaleOption = "linear"
        self.startTemperature       = 10      
        #Parameters that can be usefull to change for umbrella sampling
        self.maxItQC                = 250
        self.energyTolQC            = 2.0e-4
        self.densityTol             = 2.0e-8    
        #Setting parameters based on information that we collected on the instance construction
        self.RNG                    = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( self.seed ) )
        if not os.path.exists(_baseFolder):
            os.makedirs(_baseFolder)
    #===============================================================================    
    def ChangeDefaultParameters(self,_parameters):
        '''
        Class method to set more specifc parameters.
        
        For Bachega to think:
            there are some of those parameters that the commom user will not have interest in change
            and then you may want to not acces such method through EasyHybrid GUI.
            
        '''     
        if 'temperature'                in _parameters:
            self.temperature            =  _parameters['temperature']   
        if 'coll_freq'                  in _parameters:
            self.collFreq               = _parameters['coll_freq']
        if 'pressure_coupling'          in _parameters:
            self.pressureCoupling       = _parameters['pressure_coupling'] 
        if 'temperature_scale'          in _parameters:
            self.temperatureScaleFreq   = _parameters['temperature_scale']
        if 'timeStep'                   in _parameters:
            self.timeStep               = _parameters['timeStep']
        if 'pressureControl'            in _parameters: # the presence of the ket activates this boolean flag
            self.pressureControl        = True
        if 'sampling_Factor'             in _parameters: 
            self.samplingFactor         = _parameters['sampling_Factor']
        if 'log_frequency'              in _parameters:
            self.logFrequency           = _parameters['log_frequency']
        if 'temperature_scale_option'   in _parameters:
            self.temperatureScaleOption = _parameters['temperature_scale_option']


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

        if not os.path.exists( self.trajectoryNameCurr ):
            os.makedirs( self.trajectoryNameCurr )
        #---------------------------------------------------------------------------------------------
        VelocityVerletDynamics_SystemGeometry(self.molecule                                           ,
                                            logFrequency              = self.logFreq                  ,
                                            normalDeviateGenerator    = self.RNG                      ,
                                            steps                     = self.nsteps                   ,
                                            timeStep                  = self.timeStep                 ,
                                            trajectories              = [(trajectory,self.nsteps/100)],
                                            temperatureScaleFrequency = self.temperatureScaleFreq     ,
                                            temperatureScaleOption    = self.temperatureScaleOption   ,
                                            temperatureStart          = self.startTemperature         ,
                                            temperatureStop           = self.temperature              )
        #..............................................................................................
    
    #===================================================================================
    def RunEquilibration(self,_equiSteps):
        '''
        Run a molecular dynamics simulation for equilibration of the system
        '''

        self.nsteps             = _equiSteps
        self.trajectoryNameCurr = self.trajectoryNameEqui
    
        if not os.path.exists( self.trajectoryNameCurr ):
            os.makedirs( self.trajectoryNameCurr )

        if self.algorithm == "Verlet":
            self.runVerlet()
        elif self.algorithm == "LeapFrog":
            self.runLeapFrog()
        elif self.algorithm == "Langevin":
            self.runLangevin()
            
    #=====================================================================================    
    def RunProduction(self,_prodSteps,):
        '''
        Run a molecular dynamics simulation for data collection.
        '''
        self.nsteps             = _prodSteps
        self.trajectoryNameCurr = self.trajectoryNameProd         
        
        if not os.path.exists( self.trajectoryNameCurr ):
            os.makedirs( self.trajectoryNameCurr )

        if self.algorithm == "Verlet":
            self.runVerlet()
        elif self.algorithm == "LeapFrog":
            self.runLeapFrog()
        elif self.algorithm == "Langevin":
            self.runLangevin()

        if self.outputDCD:
            Duplicate(self.trajectoryNameCurr,self.trajectoryNameCurr+".dcd",self.molecule)
  
   
    #=====================================================================================
    def RunProductionRestricted(self,_equiSteps,_prodSteps,_samplingFactor):
        '''
        Run a simulation with the system having soft constrains defined.
        Designed for Umbrella Sampling and Steered Molecular Dynamics routines.
        '''
        self.nsteps             = _prodSteps
        self.samplingFactor     = _samplingFactor
        self.softConstraint     = True
        self.outputDCD          = False
        
        self.RunEquilibration(_equiSteps)
        self.RunProduction(_prodSteps)       
    

    #======================================================================================
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
        
        VelocityVerletDynamics_SystemGeometry(self.molecule                             ,
                                logFrequency                = self.logFreq              ,
                                normalDeviateGenerator      = self.RNG                  ,
                                steps                       = self.nsteps               ,
                                timeStep                    = self.timeStep             ,
                                temperatureScaleFrequency   = self.temperatureScaleFreq ,
                                temperatureScaleOption      = "constant"                ,
                                trajectories                = trajectory_list           ,
                                temperatureStart            = self.temperature          )

    #=====================================================================================
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
        LangevinDynamics_SystemGeometry ( self.molecule                         ,
                                          collisionFrequency     = self.collFreq,
                                          logFrequency           = self.logFreq ,
                                          normalDeviateGenerator = self.RNG     ,
                                          steps                  = self.nsteps  ,
                                          temperature            =   300.0      ,
                                          timeStep               =   0.001      ,
                                          trajectories = trajectory_list        )

    #=====================================================================================
         

#===============================================================================#
#. End of class MD =============================================================#
#===============================================================================#
       
    
    