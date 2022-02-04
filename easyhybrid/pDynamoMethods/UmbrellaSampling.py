#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#FILE = UmbrellaSampling.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#==============================================================================

#-----------------------------------------------------
import os, glob, sys
#-----------------------------------------------------
from commonFunctions import *
from MolecularDynamics import MD 
from commonFunctions import GetAtomicMass
#-----------------------------------------------------
import pymp
import numpy as np 
import matplotlib.pyplot as plt
#-----------------------------------------------------
from pMolecule import *
from pMolecule.QCModel import *
from scipy.interpolate import griddata


#********************************************************************************
class US:
    '''
    Class for setup and execute Umbrella Sampling simulations 
    ''' 
    #----------------------------------------------------------------------------   
    def __init__(self,_system,_baseFolder,_equiSteps,_prodSteps,mdMethod,RESTART=False,ADAPTATIVE=False):
        '''
        Class constructor
        '''
        self.baseName           = _baseFolder
        self.inputTraj          = " " #folder containing the pkls of the starting geometries
        self.molecule           = _system 
        self.nDim               = 0
        self.atoms              = [] # indices of the atoms 
        self.nprocs             = 1
        self.text               = " "
        self.forceC             = 600.0
        self.nsteps             = [ 1, 1 ]
        self.prodNsteps         = _prodSteps
        self.equiNsteps         = _equiSteps
        self.temperature        = 300.15
        self.maxIt              = 30
        self.rmsGT              = 0.1
        self.multipleDistance   = [ False,False ]
        self.massConstraint     = True
        self.sigma_a1_a3        = [ 0.0,0.0]
        self.sigma_a3_a1        = [ 0.0,0.0]
        self.mdMethod           = mdMethod
        self.bins               = 0
        self.restart            = RESTART

        #qc/mm parameters 
        self.maxItQC            = 1000
        self.energyTolQC        = 1e-04
        self.densityTol         = 1e-08        
        self.timeStep           = 0.001
        self.mdParameters = { }        

    #-----------------------------------------------------------------
    def ChangeDefaultParameters(self,_parameters):
    

        if 'temperature'        in _parameters:
            self.temperature    = _parameters['temperature']
        if 'maxIterations_QC'   in _parameters:
            self.maxItQC        = _parameters['maxIterations_QC']
        if 'energy_tolerance'   in _parameters:
            self.energyTolQC    = _parameters['energy_tolerance']
        if 'density_tolerancen' in _parameters:
            self.densityTol     = _parameters['density_tolerance']
        if 'timeStep'           in _parameters:
            self.timeStep       = _parameters['timeStep']
        if 'NmaxThreads'        in _parameters:
            self.nprocs         = _parameters['NmaxThreads']

        self.mdParameters       = {                                 
                                'temperature': self.temperature,\
                                'maxIterations_QC':self.maxItQC,\
                                'density_tolerancen':self.densityTol,\
                                'timeStep':self.timeStep,\
                                'energy_tolerance':self.energyTolQC
                            }

    #==========================================================================
    def SetMode(self,_atoms,_massConstraint):
        '''
        Class method to setup modes to be restrained
        '''
        #----------------------------------------------------------------------
        ndim = self.nDim # temp var to hold the index of the curren dim
        self.nDim += 1
        self.atoms.append(_atoms)
        self.sigma_a1_a3[ndim]  =  1.0 
        self.sigma_a3_a1[ndim]  = -1.0
        self.massConstraint     = _massConstraint

        #-------------------------------------------------------------------
        if len(_atoms) == 3:
            self.multipleDistance[ndim] = True

        #.---------------------------
        if self.multipleDistance[ndim]:
            #.----------------------
            if self.massConstraint:
                atomic_n1 = self.molecule.atoms.items[ self.atoms[ndim][0] ].atomicNumber
                atomic_n3 = self.molecule.atoms.items[ self.atoms[ndim][2] ].atomicNumber
                mass_a1 = GetAtomicMass(atomic_n1)
                mass_a3 = GetAtomicMass(atomic_n3)
                self.sigma_a1_a3[ndim] = mass_a1 /(mass_a1+mass_a3)
                self.sigma_a3_a1[ndim] = mass_a3 /(mass_a1+mass_a3)
                self.sigma_a3_a1[ndim] = self.sigma_a3_a1[ndim]*-1

    #============================================================================
    def SetAdaptativeParameters(self):
        pass
    
    #============================================================================
    def ChangeConvergenceParameters(self):
        pass


    #============================================================================  
    def Run1DSampling(self,_trajFolder,_sample):
        '''
        Class method to execute one-dimensional sampling
        '''
        #-----------------------------------------------
        self.inputTraj = _trajFolder
        #-----------------------------------------------
        atom1 = self.atoms[0][0]
        atom2 = self.atoms[0][1]
        atom3 = 0

        if len(self.atoms[0]) == 3:
            atom3 = self.atoms[0][2]
        #----------------------------------------------        
        weight1 = self.sigma_a3_a1[0]
        weight2 = self.sigma_a1_a3[0]        
        #---------------------------------------------
        #Adicionar outras possibilidades de carregar cordenadas
        pkl_path   = os.path.join( _trajFolder, "")
        file_lists = glob.glob( pkl_path+"*.pkl" )
        self.bins  = len(file_lists)
        #---------------------------------------------
        distance   = 0.0
        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )

        #...................................................................................
        if self.multipleDistance[0]:            
            with pymp.Parallel(self.nprocs) as p:
                for i in p.range( len(file_lists) ):
                    self.molecule.coordinates3 = ImportCoordinates3(file_lists[i])
                    #------------------------------------------------------------
                    dist_12 = self.molecule.coordinates3.Distance( atom1, atom2 )
                    dist_23 = self.molecule.coordinates3.Distance( atom2, atom3 )
                    distance = ( weight1 * dist_12) - ( weight2 * dist_23*-1)
                    #------------------------------------------------------------
                    rmodel = RestraintEnergyModel.Harmonic( distance, self.forceC )
                    restraint = RestraintMultipleDistance.WithOptions(energyModel = rmodel,  distances= [ [ atom2, atom1, weight1 ], [ atom2, atom3, weight2 ] ]) 
                    restraints['M1'] = restraint
                    #------------------------------------------------------------
                    temp = file_lists[i][:-4]
                    temp = os.path.basename(temp)
                    mdfolder = os.path.join( self.baseName, temp )
                    #------------------------------------------------------------
                    mdRun = MD(self.molecule,mdfolder,self.mdMethod)
                    mdRun.ChangeDefaultParameters(self.mdParameters)
                    mdRun.RunProductionRestricted(self.equiNsteps,self.prodNsteps,_sample)
                    #..................................................................... 
        #----------------------------------------------------------------------------------
        else:
            with pymp.Parallel(self.nprocs) as p:
                for i in p.range( len(file_lists) ):
                    self.molecule.coordinates3 = ImportCoordinates3(file_lists[i])
                    #------------------------------------------------------------                
                    distance    = self.molecule.coordinates3.Distance( atom1, atom2 )
                    rmodel      = RestraintEnergyModel.Harmonic( distance, self.forceC )
                    restraint   = RestraintDistance.WithOptions( energyModel = rmodel, point1=atom1, point2=atom2 ) 
                    restraints['M2'] = restraint
                    #------------------------------------------------------------
                    temp = file_lists[i][:-4]
                    temp = os.path.basename(temp)
                    mdfolder = os.path.join( self.baseName, temp )
                    #------------------------------------------------------------
                    mdRun = MD(self.molecule,mdfolder,self.mdMethod)
                    mdRun.ChangeDefaultParameters(self.mdParameters)
                    mdRun.RunProductionRestricted(self.equiNsteps,self.prodNsteps,_sample) 
                    #............................................................

        self.molecule.DefineRestraintModel(None) 

    #==============================================================================
    def Run2DSampling(self,_trajFolder,_sample):
        '''
        Class method to execute two-dimesninal sampling 
        ''' 
        self.samplingFactor = _sample

        pkl_path        = os.path.join( _trajFolder, "")
        self.file_lists = glob.glob( pkl_path+"*.pkl" )
        self.bins       = len(self.file_lists)
        self.mdPaths    = []

        for i in range( len(self.file_lists) ):
            coordinate_file = self.file_lists[i]
            temp    = coordinate_file[:-4]
            temp    = os.path.basename(temp)
            md_path = os.path.join(self.baseName, temp )
            self.mdPaths.append(md_path)
             
        if self.restart:   
            for fl in self.mdPaths:
                if os.path.exists( fl ):
                    self.mdPaths.remove( fl )

        #-----------------------------------------------
        if self.multipleDistance[0] and self.multipleDistance[1]:
            self.Run2DMultipleDistance()            
        elif self.multipleDistance[0] and self.multipleDistance[1] == False:            
            self.Run2DMixedDistance()
        else:
            self.Run2DSimpleDistance()  
            
    #===========================================================================================
    def Run2DMultipleDistance(self):
        '''
        '''
        atom1 = self.atoms[0][0]
        atom2 = self.atoms[0][1]
        atom3 = self.atoms[0][2]
        atom4 = self.atoms[1][0]
        atom5 = self.atoms[1][1]
        atom6 = self.atoms[1][2]

        weight1 = self.sigma_a3_a1[0]
        weight2 = self.sigma_a1_a3[0]
        weight3 = self.sigma_a3_a1[1]
        weight4 = self.sigma_a1_a3[1]

        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )


        with pymp.Parallel(self.nprocs) as p:
            for i in p.range ( self.bins) :  
                #--------------------------------------------------------
                #First confirm if the folder already exists in cases of restart
                self.molecule.coordinates3 = ImportCoordinates3( self.file_lists[i], log=None )
                #--------------------------------------------------------
                dist12      = self.molecule.coordinates3.Distance( atom1, atom2 )
                dist23      = self.molecule.coordinates3.Distance( atom2, atom3  )
                distance_1  = ( weight1 * dist12) - ( weight2 * dist23*-1)
                #--------------------------------------------------------
                rmodel      =  RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                restraint   =  RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = [ [ atom2, atom1, weight1 ],[ atom2, atom3, weight2 ] ] )
                restraints["M1"] = restraint
                #--------------------------------------------------------               
                dist45      = self.molecule.coordinates3.Distance( atom4, atom5 )
                dist56      = self.molecule.coordinates3.Distance( atom5, atom6  )
                distance_2  = ( weight1 * dist45) - ( weight2 * dist56*-1)
                #--------------------------------------------------------
                rmodel2     = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                restraint   = RestraintMultipleDistance.WithOptions( energyModel = rmodel2,  distances = [ [ atom5, atom4, weight3 ],[ atom5, atom6, weight4 ] ] )
                restraints["M2"] = restraint 
                #--------------------------------------------------------
                mdRun = MD(self.molecule,self.mdPaths[i],self.mdMethod)
                mdRun.ChangeDefaultParameters(self.mdParameters)
                mdRun.RunProductionRestricted(self.equiNsteps,self.prodNsteps,self.samplingFactor) 
        
        #.....................................................................
        self.molecule.DefineRestraintModel(None) 
        #---------------------------------------            
    
    #==========================================================================================           
    def Run2DMixedDistance(self):
        '''

        '''
        atom1 = self.atoms[0][0]
        atom2 = self.atoms[0][1]
        atom3 = self.atoms[0][2]
        atom4 = self.atoms[1][0]
        atom5 = self.atoms[1][1]

        weight1 = self.sigma_a3_a1[0]
        weight2 = self.sigma_a1_a3[0]

        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )
        
        with pymp.Parallel(self.nprocs) as p:
            for i in p.range (self.bins):                 
                #------------------------------------------------------------------------           
                self.molecule.coordinates3 = ImportCoordinates3( self.file_lists[i],log = None )
                #------------------------------------------------------------------------
                dist12      = self.molecule.coordinates3.Distance( atom1, atom2 )
                dist23      = self.molecule.coordinates3.Distance( atom2, atom3  )
                distance_1  = ( weight1 * dist12) - ( weight2 * dist23*-1)
                rmodel      =  RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                restraint   =  RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = [ [ atom2, atom1, weight1 ],[ atom2, atom3, weight2 ] ] )
                restraints["M1"] = restraint
                #-----------------------------------------------------------------------         
                distance_2  = self.molecule.coordinates3.Distance( atom4, atom5 )
                rmodel      = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                restraint   = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom4, point2= atom5)
                restraints["M2"] = restraint 
                #-----------------------------------------------------------------------
                mdRun = MD(self.molecule,self.mdPaths[i],self.mdMethod)
                mdRun.ChangeDefaultParameters(self.mdParameters)
                mdRun.RunProductionRestricted(self.equiNsteps,self.prodNsteps,self.samplingFactor)
        #---------------------------------------
        self.molecule.DefineRestraintModel(None)
        #.......................................

    #==========================================================================================
    def Run2DSimpleDistance(self):
        '''
        '''
        atom1 = self.atoms[0][0]
        atom2 = self.atoms[0][1]
        atom3 = self.atoms[1][0]
        atom4 = self.atoms[1][1]

        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )
        
        with pymp.Parallel(self.nprocs) as p:
            for i in p.range ( self.bins) :                
                #------------------------------------------------------------------------           
                self.molecule.coordinates3 = ImportCoordinates3( self.file_lists[i],log=None )
                #------------------------------------------------------------------------
                distance_1       = self.molecule.coordinates3.Distance( atom1, atom2 )                
                rmodel           = RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                restraint        = RestraintDistance.WithOptions( energyModel = rmodel, point1= atom1, point2= atom2 )
                restraints["M1"] = restraint
                #-----------------------------------------------------------------------         
                distance_2       = self.molecule.coordinates3.Distance( atom3, atom4 )
                rmodel           = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                restraint        = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom3, point2= atom4)
                restraints["M2"] = restraint  
                #-----------------------------------------------------------------------
                mdRun = MD(self.molecule,self.mdPaths[i],self.mdMethod)
                mdRun.ChangeDefaultParameters(self.mdParameters)
                mdRun.RunProductionRestricted(self.equiNsteps,self.prodNsteps,self.samplingFactor)
        #---------------------------------------
        self.molecule.DefineRestraintModel(None)
        #.......................................

    #===========================================================================================
    def Analysis(self):
        '''
        Perform analysis and plots from the production trajectories of the umbrella sampling
        '''
        pass



#==================================================================================#
#================================END OF THE CLASS==================================#
#==================================================================================#           






