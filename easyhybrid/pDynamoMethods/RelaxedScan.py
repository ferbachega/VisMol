#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#FILE = RelaxedScan.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#==============================================================================

import pymp
from commonFunctions import *
from pMolecule import *
from pMolecule.QCModel import *
from GeometrySearcher import * 

from scipy.interpolate import griddata
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.colors as colors
from matplotlib.colors import BoundaryNorm

#*****************************************************************************
class SCAN:
    '''
    Class to setup and execute relaxed surface scan procedure
    '''
    #---------------------------------------------------------------
    def __init__(self,_system,_baseFolder,_optimizer,ADAPTATIVE=False):
        '''
        Class constructor
        '''
        self.baseName           = _baseFolder
        self.molecule           = _system 
        self.nDim               = 0
        self.reactionCoordinate1= []
        self.reactionCoordinate2= []
        self.atoms              = [] 
        self.nprocs             = NmaxThreads
        self.textLog            = " "
        self.energies           = []
        self.energiesMatrix     = None
        self.DMINIMUM           = [ 0.0, 0.0 ]
        self.DINCREMENT         = [ 0.0, 0.0 ]
        self.forceC             = 2500.0
        self.massConstraint     = True
        self.multipleDistance   = [False,False]
        self.nsteps             = [ 1, 1 ]
        self.logFreq            = 50
        self.maxIt              = 800
        self.rmsGT              = 0.1
        self.optmizer           = _optimizer
        self.sigma_a1_a3        = [0.0,0.0]
        self.sigma_a3_a1        = [0.0,0.0]
        self.real_distance_1    = []
        self.real_distance_2    = []
        self.adaptative         = ADAPTATIVE
        self.text               = ""
        self.EnergyRef          = 0.0
        self.forceCRef          = self.forceC
        if not os.path.exists( os.path.join( self.baseName, "ScanTraj.ptGeo" ) ):
            os.makedirs(  os.path.join( self.baseName, "ScanTraj.ptGeo" ) )

        #set the parameters dict for the geometry search classes
        self.GeoOptPars =   { 
                            "method": self.optmizer           ,\
                            "logFrequency": self.logFreq      ,\
                            "maxIterations":self.maxIt        ,\
                            "rmsGradient": self.rmsGT
                            }
    #===========================================================================================
    def ChangeDefaultParameters(self,_parameters):
        '''
        Class method to alter deafult parameters
        '''
        #-----------------------------------------------------------
        if "rmsGradient" in _parameters:
            self.rmsGT = _parameters['rmsGradient']
        if "maxIrerations" in _parameters:
            self.maxIt = _parameters['maxIterations']
        if "log_frequency" in _parameters:
            self.logFreq = _paremeters["log_frequency"]
        if "NmaxThreads" in _parameters:
            self.nprocs = _parameters["NmaxThreads"]
        if "force_constant" in _parameters:
            self.forceC = _parameters["force_constant"]

        #-----------------------------------------------------------
        self.GeoOptPars =   { 
                            "method": self.optmizer           ,\
                            "logFrequency": self.logFreq      ,\
                            "saveTraj": "True"               ,\
                            "maxIterations":self.maxIt        ,\
                            "rmsGradient": self.rmsGT
                            }
    
    #===========================================================================================
    def ChangeConvergenceParameters(self):
        '''
        '''
        En = self.molecule.Energy()
        delta = En - self.EnergyRef

        if delta < 150.0:
            self.forceC = self.forceCRef
            self.molecule.qcModel.converger.energyTolerance  = 0.0001
            self.molecule.qcModel.converger.densityTolerance = 3e-08
            self.molecule.qcModel.converger.diisDeviation    = 1e-06
        elif delta >= 150.0:
            self.forceC = self.forceCRef - self.forceCRef*0.40
            self.molecule.qcModel.converger.energyTolerance  = 0.0003
            self.molecule.qcModel.converger.densityTolerance = 3e-08
            self.molecule.qcModel.converger.diisDeviation    = 1e-06
            if delta > 160.0 and delta < 170.0:
                self.forceC = self.forceCRef - self.forceCRef*0.50
                self.molecule.qcModel.converger.energyTolerance  = 0.0006
                self.molecule.qcModel.converger.densityTolerance = 1e-07
                self.molecule.qcModel.converger.diisDeviation    = 2e-06
            elif delta > 170.0 and delta <180.0 :
                self.forceC = self.forceCRef - self.forceCRef*0.50
                self.molecule.qcModel.converger.energyTolerance  = 0.001
                self.molecule.qcModel.converger.densityTolerance = 3e-07
                self.molecule.qcModel.converger.diisDeviation    = 5e-06
            elif delta > 180.0 and delta < 185.0:
                self.forceC = self.forceCRef - self.forceCRef*0.70
                self.molecule.qcModel.converger.energyTolerance  = 0.0015
                self.molecule.qcModel.converger.densityTolerance = 1e-06
                self.molecule.qcModel.converger.diisDeviation    = 1e-05
            elif delta > 185.0 and delta <200.0:
                self.forceC = self.forceCRef - self.forceCRef*0.70
                self.molecule.qcModel.converger.energyTolerance  = 0.003
                self.molecule.qcModel.converger.densityTolerance = 1e-05
                self.molecule.qcModel.converger.diisDeviation    = 5e-05
            elif delta > 200.0:
                self.forceC = self.forceCRef - self.forceCRef*0.70
                self.molecule.qcModel.converger.energyTolerance  = 0.003
                self.molecule.qcModel.converger.densityTolerance = 1e-04
                self.molecule.qcModel.converger.diisDeviation    = 5e-04

    #=============================================================================================
    def SetReactionCoord(self,_RC):
        '''
        Set reaction coordinate, determining initial parameters from the atoms information
        '''
        #------------------------------------------------------------
        ndim = self.nDim # temp var to hold the index of the curren dim
        self.nDim += 1
        self.atoms.append(_RC.atoms)

        self.DINCREMENT[ndim]       = _RC.increment
        self.sigma_a1_a3[ndim]      = _RC.weight13
        self.sigma_a3_a1[ndim]      = _RC.weight31
        self.DMINIMUM[ndim]         = _RC.minimumD
        self.massConstraint         = _RC.massConstraint

        if len(_RC.atoms) == 3:
            self.multipleDistance[ndim] = True
            

    #===============================================================================================
    def RunONEDimensionSCAN(self,_nsteps):
        '''
        Execute the relaxed scan with one reaction coordinate
        '''
        #-------------------------------------------------------------------------
        #Setting some local vars to ease the notation in the pDynamo methods
        #----------------------------------
        atom1 = self.atoms[0][0]
        atom2 = self.atoms[0][1]
        atom3 = 0

        if len(self.atoms[0]) == 3:
            atom3 = self.atoms[0][2]
        #----------------------------------
        weight1 = self.sigma_a1_a3[0]
        weight2 = self.sigma_a3_a1[0]              
        #---------------------------------
        self.text += "x RC1 Energy\n" 

        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )
        
        #----------------------------------------------------------------------------------------
        if self.multipleDistance[0]:
            for i in range(0,_nsteps):
                distance = self.DMINIMUM[0] + ( self.DINCREMENT[0] * float(i) ) 
                #--------------------------------------------------------------------
                rmodel = RestraintEnergyModel.Harmonic( distance, self.forceC )
                restraint = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom2, atom1, weight1 ], [ atom2, atom3, weight2 ] ] )
                restraints["RC1"] =  restraint            
                #--------------------------------------------------------------------
                relaxRun = GeometrySearcher(self.molecule, self.baseName  )
                relaxRun.ChangeDefaultParameters(self.GeoOptPars)
                relaxRun.Minimization(self.optmizer)
                #--------------------------------------------------------------------
                if i == 0:
                    en0 = self.molecule.Energy()
                    self.energies.append(0.0)
                else:
                    self.energies.append( self.molecule.Energy() - en0 )
                self.reactionCoordinate1.append( self.molecule.coordinates3.Distance( atom1 , atom2  ) - self.molecule.coordinates3.Distance( atom2, atom3  ) ) 
                Pickle( os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}.pkl".format(i) ), self.molecule.coordinates3 )
                self.text += "{} {} {} \n".format( i, self.reactionCoordinate1[i],self.energies[i]) 

        #..........................................................................................
        else:
             for i in range(_nsteps):       
                distance = self.DMINIMUM[0] + ( self.DINCREMENT[0] * float(i) )
                #--------------------------------------------------------------------
                rmodel = RestraintEnergyModel.Harmonic( distance, self.forceC )
                restraint = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom1, point2= atom2 )
                restraints["RC1"] =  restraint            
                #--------------------------------------------------------------------
                relaxRun = GeometrySearcher(self.molecule,self.baseName)
                relaxRun.ChangeDefaultParameters(self.GeoOptPars)
                relaxRun.Minimization(self.optmizer)
                #--------------------------------------------------------------------
                if i == 0:
                    en0 = self.molecule.Energy()
                    self.energies.append(0.0)
                else:
                    self.energies.append( self.molecule.Energy() - en0 )
                #--------------------------------------------------------------------
                self.reactionCoordinate1.append( self.molecule.coordinates3.Distance( atom1 , atom2  ) )  
                Pickle( os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}.pkl".format(i) ), self.molecule.coordinates3 ) 
                self.text += "{} {} {} \n".format( i, self.reactionCoordinate1[i], self.energies[i] )
        #---------------------------------------
        self.molecule.DefineRestraintModel(None)
    #===================================================================================================
    def Run2DScan(self,_nsteps_x,_nsteps_y):
        '''
        '''
        #------------------------------------------------------
        self.text += "x y RC1 RC2 Energy\n" 
        #------------------------------------------------------
        restraints = RestraintModel( )
        self.molecule.DefineRestraintModel( restraints )
        #------------------------------------------------------
        distance_1 = 0.0 
        distance_2 = 0.0
        #------------------------------------------------------
        X = _nsteps_x
        Y = _nsteps_y
        self.nsteps[0] = X       
        self.nsteps[1] = Y  
        self.energiesMatrix = pymp.shared.array( (X,Y), dtype=float ) 
        self.reactionCoordinate1 = pymp.shared.array( (X,Y), dtype=float )   
        self.reactionCoordinate2 = pymp.shared.array( (X,Y), dtype=float )   
        #-----------------------------------------------------------------------------------------------
        #Define the origin point of the relaxed surface scan, AKA the 0,0 point
        coordinateFile = os.path.join( self.baseName ,"ScanTraj.ptGeo","frame{}_{}.pkl".format( 0, 0 ) )
        relaxRun = GeometrySearcher( self.molecule, self.baseName )
        relaxRun.ChangeDefaultParameters( self.GeoOptPars )
        relaxRun.Minimization(self.optmizer)
        #-----------------------------------------------------------------------------------------------
        self.EnergyRef = self.en0 = self.molecule.Energy(log=None)         
        Pickle( coordinateFile, self.molecule.coordinates3 )


        if self.multipleDistance[0] and self.multipleDistance[1]:
            self.Run2DScanMultipleDistance(X,Y)            
        elif self.multipleDistance[0] and self.multipleDistance[1] == False:            
            self.Run2DMixedDistance(X,Y)
        else:
            self.Run2DSimpleDistance(X,Y)

        for i in range(X):
            for j in range(Y):
                self.text += "{} {} {} {} {}\n".format( i,j,self.reactionCoordinate1[i,j], self.reactionCoordinate2[i,j], self.energiesMatrix[i,j])
    #=======================================================
    def Run2DSimpleDistance(self, X, Y ):
        '''
        '''
        atom1 = self.atoms[0][0]
        atom2 = self.atoms[0][1]
        atom3 = self.atoms[1][0]
        atom4 = self.atoms[1][1]

        restraints = RestraintModel( )
        self.molecule.DefineRestraintModel( restraints )

        self.reactionCoordinate1[ 0,0 ] = self.molecule.coordinates3.Distance( atom1, atom2 ) 
        self.reactionCoordinate2[ 0,0 ] = self.molecule.coordinates3.Distance( atom3, atom4 ) 

        with pymp.Parallel(self.nprocs) as p:
            for i in p.range ( 1, X ):  
                #.---------------------------------------------------------------------------------------------             
                distance_1 = self.DMINIMUM[0] + ( self.DINCREMENT[0] * float(i) ) 
                rmodel     =  RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                restraint  =  RestraintDistance.WithOptions( energyModel = rmodel,  point1=atom1, point2=atom2  )
                restraints["RC1"] = restraint                
                #----------------------------------------------------------------------------------------------                
                distance_2  = self.DMINIMUM[1]
                rmodel      = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                restraint   = RestraintDistance.WithOptions( energyModel = rmodel, point1=atom3, point2=atom4 )                    
                restraints["RC2"] = restraint  
                #----------------------------------------------------------------------------------------------                   
                initCoordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( 0 , 0) ) 
                self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile, log=None )
                if self.adaptative:
                    self.ChangeConvergenceParameters()
                #----------------------------------------------------------------------------------------------
                coordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( i, 0 ) )
                relaxRun = GeometrySearcher( self.molecule, self.baseName )
                relaxRun.ChangeDefaultParameters(self.GeoOptPars)
                relaxRun.Minimization(self.optmizer)
                #----------------------------------------------------------------------------------------------
                self.energiesMatrix[ i,0 ]      = self.molecule.Energy(log=None) - self.en0
                self.reactionCoordinate1[ i,0 ] = self.molecule.coordinates3.Distance( atom1, atom2 ) 
                self.reactionCoordinate2[ i,0 ] = self.molecule.coordinates3.Distance( atom3, atom4 )   
                #-----------------------------------------------------------------------------------
                Pickle( coordinateFile, self.molecule.coordinates3 )
                #....................................................
            #-------------------------------------------------------------------------------------------
        with pymp.Parallel(self.nprocs) as p:
            #Pergomr the calculations for the rest of the grid
            for i in p.range ( 0, X ):
                #----------------------------------------------------------------------------------------------
                distance_1 = self.DMINIMUM[0] + ( self.DINCREMENT[0] * float(i) )
                rmodel     = RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                restraint  = RestraintDistance.WithOptions(energyModel =rmodel, point1=atom1, point2=atom2  )
                restraints["RC1"] = restraint
                #----------------------------------------------------------------------------------------------
                for j in range( 1, Y ):

                    distance_2  = self.DMINIMUM[1] + ( self.DINCREMENT[1] * float(j) )
                    rmodel      = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                    restraint   = RestraintDistance.WithOptions(energyModel = rmodel, point1=atom3, point2=atom4  )
                    restraints["RC2"] = restraint                    
                    #------------------------------------------------------------------------------------------
                    initCoordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo" , "frame{}_{}.pkl".format( i, j-1 ) )                    
                    self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile, log=None )             
                    if self.adaptative:
                        self.ChangeConvergenceParameters()
                    #----------------------------------------------------------------------------------------------
                    coordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( i, j ) )
                    relaxRun = GeometrySearcher( self.molecule, self.baseName  )
                    relaxRun.ChangeDefaultParameters( self.GeoOptPars )
                    relaxRun.Minimization(self.optmizer) 
                    #----------------------------------------------------------------------------------------------
                    self.energiesMatrix[ i,j ] = self.molecule.Energy(log=None) - self.en0
                    self.reactionCoordinate1[ i,j ] = self.molecule.coordinates3.Distance( atom1, atom2 ) 
                    self.reactionCoordinate2[ i,j ] = self.molecule.coordinates3.Distance( atom3, atom4 )
                    #-----------------------------------------------------------------------------------
                    Pickle( coordinateFile, self.molecule.coordinates3 )
                    #-----------------------------------------------------------------------------------
        self.molecule.DefineRestraintModel(None)
        #-----------------------------------------------------------------------------------------------

    #=======================================================   
    def Run2DMixedDistance(self,X, Y ):
        '''
        '''
        atom1 = self.atoms[0][0]
        atom2 = self.atoms[0][1]
        atom3 = self.atoms[0][2]
        atom4 = self.atoms[1][0]
        atom5 = self.atoms[1][1]
        
        weight1 = self.sigma_a1_a3[0]
        weight2 = self.sigma_a3_a1[0]

        restraints = RestraintModel( )
        self.molecule.DefineRestraintModel( restraints )
        
        self.reactionCoordinate1[ 0,0 ] = self.molecule.coordinates3.Distance( atom1, atom2 ) - self.molecule.coordinates3.Distance( atom3, atom2 )
        self.reactionCoordinate2[ 0,0 ] = self.molecule.coordinates3.Distance( atom4, atom5 )
       
        with pymp.Parallel(self.nprocs) as p:
            for i in p.range ( 1, X ):  
                #--------------------------------------------------------------------------------  
                distance_1 = self.DMINIMUM[0] + self.DINCREMENT[0] * float(i)
                rmodel     = RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                restraint  = RestraintMultipleDistance.WithOptions(energyModel = rmodel, distances = [ [ atom2, atom1, weight1 ],[ atom2, atom3, weight2 ] ] )
                restraints["RC1"] = restraint                
                #--------------------------------------------------------------------------------
                distance_2  = self.DMINIMUM[1]
                rmodel      = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                restraint   = RestraintDistance.WithOptions( energyModel = rmodel, point1=atom4, point2=atom5 )                
                restraints["RC2"] = restraint  
                #---------------------------------------------------------------------------------                    
                initCoordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format(0,0) ) 
                self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile, log = None )
                #--------------------------------------------------------------------------------
                if self.adaptative:
                    self.ChangeConvergenceParameters() 
                #--------------------------------------------------------------------------------
                coordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( i, 0 ) )
                relaxRun = GeometrySearcher( self.molecule, self.baseName )
                relaxRun.ChangeDefaultParameters( self.GeoOptPars )
                relaxRun.Minimization(self.optmizer)
                #-------------------------------------------------------------------------------- 
                self.energiesMatrix[ i,0 ]      = self.molecule.Energy(log=None) - self.en0
                self.reactionCoordinate1[ i,0 ] = self.molecule.coordinates3.Distance( atom1, atom2 ) - self.molecule.coordinates3.Distance( atom3, atom2 )
                self.reactionCoordinate2[ i,0 ] = self.molecule.coordinates3.Distance( atom4, atom5 )              
                #-----------------------------------------------------------------------------------
                Pickle( coordinateFile, self.molecule.coordinates3 ) 
                #..................................................
        #---------------------------------------------------------------------------------------------
        with pymp.Parallel(self.nprocs) as p:
            #Pergomr the calculations for the rest of the grid
            for i in p.range ( 0, X ):
                distance_1 = self.DMINIMUM[0] + self.DINCREMENT[0] * float(i) 
                rmodel     =  RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                restraint  =  RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom2, atom1, weight1 ],[ atom2, atom3, weight2 ] ] )
                restraints["RC1"] = restraint
                #-----------------------------------------------------------------------------------
                
                for j in range( 1, Y ):
                    distance_2  = self.DMINIMUM[1] + ( self.DINCREMENT[1] * float(j) )
                    rmodel      = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                    restraint   = RestraintDistance.WithOptions( energyModel = rmodel, point1=atom4, point2=atom5 )
                    restraints["RC2"] = restraint  
                    #-----------------------------------------------------------------------------------
                    initCoordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo" , "frame{}_{}.pkl".format( i, j-1 ) )                 
                    #-----------------------------------------------------------------------------------
                    self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile, log=None )             
                    #-----------------------------------------------------------------------------------
                    if self.adaptative:
                        self.ChangeConvergenceParameters()
                    #-----------------------------------------------------------------------------------
                    coordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( i, j ) )
                    relaxRun = GeometrySearcher( self.molecule, self.baseName  )
                    relaxRun.ChangeDefaultParameters( self.GeoOptPars )
                    relaxRun.Minimization(self.optmizer)
                    #-----------------------------------------------------------------------------------
                    self.energiesMatrix[ i,j ]      = self.molecule.Energy(log=None) - self.en0
                    self.reactionCoordinate1[ i,j ] = self.molecule.coordinates3.Distance( atom1, atom2 ) - self.molecule.coordinates3.Distance( atom3, atom2 )
                    self.reactionCoordinate2[ i,j ] = self.molecule.coordinates3.Distance( atom4, atom5 )
                    #-----------------------------------------------------------------------------------
                    Pickle( coordinateFile, self.molecule.coordinates3 )
                    #...................................................
        self.molecule.DefineRestraintModel(None)
    
    #===========================================================
    def Run2DScanMultipleDistance(self, X, Y ):
        '''
        '''
        atom1 = self.atoms[0][0]
        atom2 = self.atoms[0][1]
        atom3 = self.atoms[0][2]
        atom4 = self.atoms[1][0]
        atom5 = self.atoms[1][1]
        atom6 = self.atoms[1][2]

        weight1 = self.sigma_a1_a3[0]
        weight2 = self.sigma_a3_a1[0]
        weight3 = self.sigma_a1_a3[1]
        weight4 = self.sigma_a3_a1[1]

        restraints = RestraintModel( )
        self.molecule.DefineRestraintModel( restraints )

        self.reactionCoordinate1[ 0,0 ] = self.molecule.coordinates3.Distance( atom1, atom2 ) - self.molecule.coordinates3.Distance( atom3, atom2 )
        self.reactionCoordinate2[ 0,0 ] = self.molecule.coordinates3.Distance( atom4, atom5 ) - self.molecule.coordinates3.Distance( atom6, atom5 )
                
        #-------------------------------------------------------------------------------------
        with pymp.Parallel(self.nprocs) as p:
            for i in p.range ( 1, X ):  
            #.---- ----------------------------------------------------------------------------            
                distance_1 = self.DMINIMUM[0] + self.DINCREMENT[0] * float(i) 
                rmodel     = RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                restraint  = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = [ [ atom2, atom1, weight1 ] , [ atom2, atom3, weight2 ] ] )
                restraints["RC1"] = restraint
                #---- ----------------------------------------------------------------------------        
                distance_2  = self.DMINIMUM[1]
                rmodel      = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                restraint   = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = [ [ atom5, atom4, weight3 ],[ atom5, atom6, weight4 ] ] )
                restraints["RC2"] = restraint  
                #---------------------------------------------------------------------------------
                initCoordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format(0,0) ) 
                self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile, log = None )             
                #-----------------------------------------------------------------------------------
                if self.adaptative:
                    self.ChangeConvergenceParameters()
                #-----------------------------------------------------------------------------------
                relaxRun = GeometrySearcher( self.molecule, self.baseName )
                relaxRun.ChangeDefaultParameters(self.GeoOptPars)
                relaxRun.Minimization(self.optmizer)
                #-----------------------------------------------------------------------------------
                self.energiesMatrix[ i,0 ]      = self.molecule.Energy(log=None) - self.en0
                self.reactionCoordinate1[ i,0 ] = self.molecule.coordinates3.Distance( atom1, atom2 ) - self.molecule.coordinates3.Distance( atom3, atom2 )
                self.reactionCoordinate2[ i,0 ] = self.molecule.coordinates3.Distance( atom4, atom5 ) - self.molecule.coordinates3.Distance( atom6, atom5 )
                #-----------------------------------------------------------------------------------
                coordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( i, 0 ) )                       
                Pickle( coordinateFile, self.molecule.coordinates3 )        
        #........................................................................................
        with pymp.Parallel(self.nprocs) as p:
            for i in p.range ( 0, X ):
                distance_1  = self.DMINIMUM[0] + self.DINCREMENT[0] * float(i) 
                rmodel      = RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                restraint   = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom2, atom1, weight1 ],[ atom2, atom3, weight2 ] ] )
                restraints["RC1"] = restraint                       
                #---------------------------------------------------------------------------------
                for j in range( 1, Y ):
                    distance_2  =  self.DMINIMUM[1] + self.DINCREMENT[1] * float(j) 
                    rmodel      = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                    restraint   = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = [ [ atom5, atom4, weight3 ],[ atom5, atom6, weight4 ] ] )
                    restraints["RC2"] = restraint  
                    #---------------------------------------------------------------------------------
                    initCoordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo" , "frame{}_{}.pkl".format( i, j-1 ) )  
                    #----------------------------------------------------------------------------              
                    self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile, log = None )             
                    if self.adaptative:
                        self.ChangeConvergenceParameters()
                    #----------------------------------------------------------------------------
                    relaxRun = GeometrySearcher( self.molecule, self.baseName  )
                    relaxRun.ChangeDefaultParameters(self.GeoOptPars)
                    relaxRun.Minimization( self.optmizer )
                    #----------------------------------------------------------------------------
                    self.energiesMatrix[ i,j ]      = self.molecule.Energy(log=None) - self.en0
                    self.reactionCoordinate1[ i,j ] = self.molecule.coordinates3.Distance( atom1, atom2 ) - self.molecule.coordinates3.Distance( atom3, atom2 )
                    self.reactionCoordinate2[ i,j ] = self.molecule.coordinates3.Distance( atom4, atom5 ) - self.molecule.coordinates3.Distance( atom6, atom5 )
                    #-----------------------------------------------------------------------------------
                    coordinateFile = os.path.join( self.baseName, "ScanTraj.ptGeo", "frame"+str(i)+"_"+str(j)+".pkl" )
                    Pickle( coordinateFile, self.molecule.coordinates3 )                    
        #--------------------------------------                
        self.molecule.DefineRestraintModel(None)

    #=======================================================================================
    def Finalize(self):
        '''
        Writing logs, making plots and saving trajectories
        '''
        #-----------------------------------------------------------------
        if self.nDim == 1:
            #..................................................
            trajNameDCD = self.baseName + ".dcd"
            trajName = os.path.join( self.baseName, "ScanTraj.ptGeo" )
            Duplicate( trajName, trajNameDCD, self.molecule )
            #..................................................
        textLog = open( self.baseName+"_SCAN{}D.log".format(self.nDim), "w" ) 
        textLog.write(self.text)
        textLog.close() 
            #..................................................
           
       

#==============================================================================#
#=====================END OF CLASS FILE========================================#
#==============================================================================#
