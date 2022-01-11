#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = GeometrySearcher.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#==============================================================================

import os, sys

#importing our library functions
import commonFunctions
from LogFile import LogFile

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

class GeometrySearcher:
    '''
    Class to handle with pDynamo methods that search geometries for the system, such as global/local minimuns
    as saddle points and reaction path trajectories. 
    '''
    #.-------------------------------------------------------------------------
    
    def __init__(self,_system,_baseFolder):
        '''
        Class constructor.
        '''
        self.molecule       = _system
        self.baseName       = _baseFolder
        self.optAlg         = "ConjugatedGradient"
        self.InitCrd3D      = _system.coordinates3
        self.finalCrd3D     = None
        self.massWeighting  = False
        self.logFreq        = 50 # deafult value for otimizations, must to be changed through the specific class method
        self.saveTraj       = False # optimization trj are not generally usefull and generate a lot of data 
        self.trajectotyName = os.path.join( _baseFolder, "trj" ) 
        self.logname        = _baseFolder + ".log"
        self.savePdb        = True
        self.traj           = None
        self.outputDCD      = True

        #Control and tolerance parameters
        self.maxIt          = 500
        self.maxItQC        = 250
        self.rmsGrad        = 0.1

    #-----------------------------------------------------------------------
    def ChengeDefaultParameters(self,_parameters):
        '''
        Class method to modify default parameters for the minimization runs
        '''
        if "log_frequency"      in _parameters:
            self.logFreq        = _paremeters["log_frequency"]
        if "not_save_pdb"       in _parameters:
            self.savePdb        = False
        if "save_traj"          in _parameters:
            self.saveTraj       = True
        if "not_save_dcd"       in _parameters:
            self.outputDCD      = False
        if "maxIrerations"      in _parameters:
            self.maxIt          = _parameters['maxIterations']
        if "maxIterations_QC"   in _parameters:
            self.maxItQC        = _parameters['maxIterations_QC']
        if 'rmsGradient'        in _parameters:
            self.rmsGrad        = _parameters['rmsGradient']

    #.-----------------------------------------------------------------------
    # Main minimization class method
    def Minimization(self,_optimizer):
        '''
        Class method to execute the minimization routine for search of geometry corresponding to local minima
        '''
        self.optAlg = _optimizer
        self.traj = ExportTrajectory( self.trajectoryName, self.molecule ) 
        
        # run the minimization for the chosen algorithm
        if self.optAlg == "ConjugatedGradient":
            self.RunConjugatedGrad()
        elif self.optAlg == "SteepestDescent":
            self.RunSteepestDescent()
        elif self.optAlg == "LFBGS":
            self.RunLFBGS()
        elif self.optAlg == "QuasiNewton":
            self.RunQuasiNewton()

        self.finalCrd3D = self.molecule.coordinates3
       
        #Save structures and/or trajectories
        if self.savePdb:
            pdbFile =  self.baseName + "_opt.pdb"


    #.-----------------------------------------------------------------------
    #Minimizers methods
    def RunConjugatedGrad(self):
        '''
        Class method to apply the conjugated gradient minimizer
        '''
        Log = TextLogFileWriter(self.logname)

        ConjugateGradientMinimize_SystemGeometry(self.molecule                          ,                
                                                 log                    = Log           ,
                                                 logFrequency           = self.logFreq  ,
                                                 trajectory             = self.traj     ,
                                                 maximumIterations      = self.maxIt    ,
                                                 rmsGradientTolerance   = self.rmsGrad  )

    #.------------------------------------------------------------------------
    def RunSteepestDescent(self):
        '''
        Class method to apply the steepest descent minimizer
        '''
        Log = TextLogFileWriter(self.logname)

        SteepestDescentMinimize_SystemGeometry(self.molecule                           ,               
                                                log                     = Log          ,
                                                logFrequency            = self.logFreq ,
                                                trajectory              = self.traj    ,
                                                maximumIterations       = self.maxIt   ,
                                                rmsGradientTolerance    = self.rmsGrad )


    #.------------------------------------------------------------------------
    def RunLFBGS(self):
        '''
        Class method to apply the LFBGS minimizer
        '''
        Log = TextLogFileWriter(self.logname)

        LBFGSMinimize_SystemGeometry(self.molecule                              ,                
                                    log  = Log                                  ,
                                    logFrequency         = self.logFreq         ,
                                    trajectory           = self.traj            ,
                                    maximumIterations    = self.maxIt           ,
                                    rmsGradientTolerance = self.rmsGrad         )
    
    #.------------------------------------------------------------------------
    def RunQuasiNewton(self):
        '''
        Class method to apply the Quaisi-Newton minimizer
        '''
        Log = TextLogFileWriter(self.logname)

        QuasiNewtonMinimize_SystemGeometry( system                              ,                
                                            log                  = Log          ,
                                            logFrequency         = self.logFreq ,
                                            trajectories         = self.traj    ,
                                            maximumIterations    = self.maxIt   ,
                                            rmsGradientTolerance = self.rmsGrad )

    #.--------------------------------------------
    # Reaction path searchers
    def NudgedElasticBand(self,_initCoord,_finalCoord,_nbins,_rmsGIS):
        '''
        Nudget Elastic Band procedure to estimate a reaction path
        '''

        _rmdGIS = _parameters["RMS_growing_intial_string"]

        self.trajectoryName = os.path.join(self.baseName + "NEB.trj")

        #Note: is interesting to think in a window were the user select the initial and final coords
        # here we excpect to ibe in pkl probably from a scan or optimization already done using the software
        self.InitCrd3D  = Unpickle( _initCoord ) # we excpect to ibe in pkl probably from a scan or optimization already done using the software
        self.finalCrd3D = Unpickle( _finalCoord ) 

        GrowingStringInitialPath (self.system ,_nBins, self.InitCrd3D, self.finalCrd3D, self.trajectoryName ,rmsGradientTolerance=_rmdGIS )

        self.traj = ExportTrajectory( self.trajectoryName, self.molecule, append=True ) 

        ChainOfStatesOptimizePath_SystemGeometry (  self.system                 , 
                                                    self.traj                   ,
                                                    logFrequency         = 1    ,
                                                    maximumIterations    = 1000 ,
                                                    fixedTerminalImages  = True ,
                                                    rmsGradientTolerance = 0.1  )
        if self.outputDCD:
            DCDTrajectory_FromSystemGeometryTrajectory( self.trajectoryName+".dcd" , self.trajectoryName, self.system )


    #.--------------------------------------------
    def SelfAvoidWalking(self,_parameters):
        '''
        Self-Avoid-Walking procedure to estimate a reaction path
        '''
        self.traj = ExportTrajectory( self.trajectoryName, self.molecule ) 


    #.--------------------------------------------
    def BakerSaddle(self,_parameters):
        '''
        Class method to search saddle-points transition structure
        '''

    #.--------------------------------------------
    #Helpers
    def CleanDir(self):
        '''
        Class method to clean unrequired files from the directory
        '''
        
    #.--------------------------------------------
    def UnitTest(self):
        '''
        '''
        pass      
  
        