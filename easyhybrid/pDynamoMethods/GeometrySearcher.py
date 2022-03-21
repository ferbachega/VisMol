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
from pScientific               import *                                     
from pScientific.Arrays        import *                                     
from pScientific.Geometry3     import *                 
from pSimulation               import *
#*********************************************************************************
class GeometrySearcher:
    '''
    Class to handle with pDynamo methods that search geometries for the system, such as global/local minimuns
    as saddle points and reaction path trajectories. 
    '''
    #.-------------------------------------------------------------------------   
    def __init__(self,_system,_baseFolder,_trajName=None):
        '''
        Class constructor.
        '''
        self.molecule       = _system
        self.baseName       = _baseFolder
        self.optAlg         = "ConjugatedGradient"
        self.InitCrd3D      = Clone(_system.coordinates3)
        self.finalCrd3D     = None
        self.massWeighting  = False
        self.logFreq        = 50 # deafult value for otimizations, must to be changed through the specific class method
        self.trajectoryName = _trajName 
        self.savePdb        = False
        self.traj           = None
        self.saveFormat     = True       
        self.rmsGrad        = 0.1
        self.maxIt          = 500
        self.saveFrequency  = 0
    #=========================================================================
    def ChangeDefaultParameters(self,_parameters):
        '''
        Class method to modify default parameters for the minimization runs
        '''       
        if "save_pdb" in _parameters:
            if _parameters["save_pdb"] == True:
                self.savePdb = True
        if "maxIterations" in _parameters:
            self.maxIt = _parameters['maxIterations']            
        if "log_frequency" in _parameters:
            self.logFreq = _parameters["log_frequency"] 
        if "save_format" in _parameters:
            self.saveFormat = _parameters["save_format"]
        if "save_frequency" in _parameters:
            self.saveFrequency = _parameters["save_frequency"]        
        if "rmsGradient" in _parameters:
            self.rmsGrad = _parameters['rmsGradient']
    #======================================================================================
    # Main minimization class method
    def Minimization(self,_optimizer):
        '''
        Execute the minimization routine for search of geometry corresponding to local minima
        '''
        #------------------------------------------------------------------
        self.optAlg = _optimizer                 
        # run the minimization for the chosen algorithm
        if   self.optAlg == "ConjugatedGradient": self.RunConjugatedGrad()
        elif self.optAlg == "SteepestDescent"   : self.RunSteepestDescent()
        elif self.optAlg == "LFBGS"             : self.RunLFBGS()
        elif self.optAlg == "QuasiNewton"       : self.RunQuasiNewton()
        elif self.optAlg == "FIRE"              : self.RunFIREmin()
        self.finalCrd3D = Clone(self.molecule.coordinates3)
    #=============================================================================
    #Minimizers methods
    def RunConjugatedGrad(self):
        '''
        Class method to apply the conjugated gradient minimizer
        '''
        if self.trajectoryName == None:
            ConjugateGradientMinimize_SystemGeometry(self.molecule                      ,                
                                                 logFrequency           = self.logFreq  ,
                                                 maximumIterations      = self.maxIt    ,
                                                 rmsGradientTolerance   = self.rmsGrad  )
        else:            
            trajectory = ExportTrajectory( self.trajectoryName, self.molecule, log=None )
            ConjugateGradientMinimize_SystemGeometry(self.molecule                        ,                
                                                 logFrequency           = self.logFreq    ,
                                                 trajectories           = [(trajectory, self.saveFrequency)],
                                                 maximumIterations      = self.maxIt      ,
                                                 rmsGradientTolerance   = self.rmsGrad    )

    #=====================================================================================
    def RunSteepestDescent(self):
        '''
        Class method to apply the steepest descent minimizer
        '''        
        if self.trajectoryName == None:
            SteepestDescentMinimize_SystemGeometry(self.molecule                       ,               
                                                logFrequency            = self.logFreq ,
                                                maximumIterations       = self.maxIt   ,
                                                rmsGradientTolerance    = self.rmsGrad )
        else:
            trajectory = ExportTrajectory( self.trajectoryName, self.molecule, log=None  )
            SteepestDescentMinimize_SystemGeometry(self.molecule                       ,               
                                                logFrequency            = self.logFreq ,
                                                trajectories            = [(trajectory, self.saveFrequency)],
                                                maximumIterations       = self.maxIt   ,
                                                rmsGradientTolerance    = self.rmsGrad )
    #============================================================================
    def RunLFBGS(self):
        '''
        Class method to apply the LFBGS minimizer
        '''        
        if self.trajectoryName == None:
            LBFGSMinimize_SystemGeometry(self.molecule                          ,                
                                    logFrequency         = self.logFreq         ,
                                    maximumIterations    = self.maxIt           ,
                                    rmsGradientTolerance = self.rmsGrad         )
        else:
            trajectory = ExportTrajectory( self.trajectoryName, self.molecule, log=None )
            LBFGSMinimize_SystemGeometry(self.molecule                          ,                
                                    logFrequency         = self.logFreq         ,
                                    trajectories         = [(trajectory, self.saveFrequency)],
                                    maximumIterations    = self.maxIt           ,
                                    rmsGradientTolerance = self.rmsGrad         )    
    #=============================================================================
    def RunQuasiNewton(self):
        '''
        Class method to apply the Quaisi-Newton minimizer
        '''        
        if self.trajectoryName == None: 
            QuasiNewtonMinimize_SystemGeometry( self.molecule                       ,                
                                                logFrequency         = self.logFreq ,
                                                maximumIterations    = self.maxIt   ,
                                                rmsGradientTolerance = self.rmsGrad )
        else:
            trajectory = ExportTrajectory( self.trajectoryName, self.molecule, log=None )
            QuasiNewtonMinimize_SystemGeometry( self.molecule                       ,                
                                                logFrequency         = self.logFreq ,
                                                trajectories         = [(trajectory, self.saveFrequency)],
                                                maximumIterations    = self.maxIt   ,
                                                rmsGradientTolerance = self.rmsGrad )
    #==============================================================================
    def RunFIREmin(self):
        '''

        '''
        if self.traj == None:
            FIREMinimize_SystemGeometry( self.molecule                  ,                
                                    logFrequency         = self.logFreq ,
                                    maximumIterations    = self.maxIt   ,
                                    rmsGradientTolerance = self.rmsGrad )
        else:
            trajectory = ExportTrajectory( self.trajectoryName, self.molecule, log=None )
            FIREMinimize_SystemGeometry( self.molecule                                       ,                
                                    logFrequency         = self.logFreq                      ,
                                    trajectories         = [(trajectory, self.saveFrequency)],
                                    maximumIterations    = self.maxIt                        ,
                                    rmsGradientTolerance = self.rmsGrad                      )        
    #=============================================================================
    # Reaction path searchers
    def NudgedElasticBand(self,_parameters):
        '''
        Nudget Elastic Band procedure to estimate a reaction path
        '''
        #-------------------------------------------------------------------------
        _rmdGIS = _parameters["RMS_growing_intial_string"]
        self.trajectoryName = os.path.join(self.baseName + "NEB.ptGeo")
        #Note: is interesting to think in a window were the user select the initial and final coords
        # here we excpect to ibe in pkl probably from a scan or optimization already done using the software
        self.InitCrd3D  = ImportCoordinates3( _parameters["init_coord"], log=None  ) # we excpect to ibe in pkl probably from a scan or optimization already done using the software
        self.finalCrd3D = ImportCoordinates3( _parameters["final_coord"], log=None ) 

        #---------------------------------------------------------------------------------
        GrowingStringInitialPath (self.molecule ,_parameters["NEB_nbins"], self.InitCrd3D, self.finalCrd3D, self.trajectoryName ,rmsGradientTolerance=_parameters["RMS_growing_intial_string"] )

        trajectory = ExportTrajectory( self.trajectoryName, self.molecule, append=True ) 

        ChainOfStatesOptimizePath_SystemGeometry (  self.molecule               , 
                                                    trajectory                  ,
                                                    logFrequency         = 1    ,
                                                    maximumIterations    = 1000 ,
                                                    fixedTerminalImages  = True ,
                                                    rmsGradientTolerance = 0.1  )        
    #========================================================================================
    def SelfAvoidWalking(self,_parameters):
        '''
        Self-Avoid-Walking procedure to estimate a reaction path
        '''
        self.traj = ExportTrajectory( self.trajectoryName, self.molecule ) 
    #========================================================================================
    def BakerSaddle(self,_parameters):
        '''
        Class method to search saddle-points transition structure
        '''
    #=========================================================================================
    def CalculateRMS(self):
        '''
        Calculate the root mean square of deviation of the final coordinate found with the first set given.
        '''
        masses = Array.FromIterable ( [ atom.mass for atom in self.molecule.atoms ] )
        self.InitCrd3D.Superimpose ( self.finalCrd3D, weights = masses )
        rms = self.InitCrd3D.RootMeanSquareDeviation ( self.finalCrd3D, weights = masses )
        print("Root Mean Sqaute of Deviation of the optimized structure from the initial: {}".format(rms))
    #===========================================================================================
    def Finalize(self):
        '''
        Finaluze the Geometry searcher procedures, save structures and/or trajectories
        '''
        self.CalculateRMS()
        #----------------------------------------------------------------------
        #Save structures and/or trajectories
        if self.savePdb:
            pdbFile = self.baseName + "opt_{}.pdb".format(self.optAlg)
            i = 0;
            while os.path.exists(pdbFile):
                pdbFile = self.baseName + "_#{}_opt_{}.pdb".format(i,self.optAlg)
                i += 1
            ExportSystem(pdbFile,self.molecule)
        #----------------------------------------------------------------------
        if not self.saveFormat == None:
            if self.saveFormat != self.trajectoryName:
                traj_save = os.path.splitext(self.trajectoryName)[0] + self.saveFormat
                print(traj_save,self.trajectoryName,self.saveFormat)
                Duplicate(self.trajectoryName,traj_save,self.molecule)
#================================================================================================#
#======================================END OF THE FILE===========================================#
#================================================================================================#
