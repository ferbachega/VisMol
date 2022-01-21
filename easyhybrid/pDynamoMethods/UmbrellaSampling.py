#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#FILE = UmbrellaSampling.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#==============================================================================

import os
from MolecularDynamics import MD 
from commonFunctions import GetAtomicMass
import glob

import pymp
from commonFunctions import *
from pMolecule import *
from pMolecule.QCModel import *

from scipy.interpolate import griddata
import numpy as np 
import matplotlib.pyplot as plt


class US:
    '''
    Class for setup and execute Umbrella Sampling simulations 
    '''    
    def __init__(self,_system,_baseFolder,_equiSteps,_prodSteps,mdMethod):
        '''
        Class constructor
        '''
        self.baseName           = _baseFolder
        self.inputTraj          = " " #folder containing the pkls of the starting geometries
        self.molecule           = _system 
        self.nDim               = 0
        self.atoms              = [] # indices of the atoms 
        self.nprocs             = 8
        self.text               = " "
        self.forceC             = 500.0
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

        self.mdParameters       = {                                 
                                'temperature': self.temperature,\
                                'maxIterations_QC':self.maxItQC,\
                                'density_tolerancen':self.densityTol,\
                                'timeStep':self.timeStep,\
                                'energy_tolerance':self.energyTolQC
                            }

    #----------------------------------------------------------------------------------
    def SetMode(self,_atoms,_massConstraint):
        '''
        Class method to setup modes to be restrained
        '''
        ndim = self.nDim # temp var to hold the index of the curren dim
        self.nDim += 1
        self.atoms.append(_atoms)
        self.sigma_a1_a3[ndim]  = 1 
        self.sigma_a3_a1[ndim]  = -1
        self.massConstraint     = _massConstraint

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

                
    #----------------------------------------------------------------------------------  
    def RunONEDimensional(self,_trajFolder,_sample):
        '''
        Class method to execute one-dimensional sampling
        '''
        atom1 = self.atoms[0][0]
        atom2 = self.atoms[0][1]
        atom3 = 0

        weight1 = self.sigma_a3_a1[0]
        weight2 = self.sigma_a1_a3[0]

        if len(self.atoms[0]) == 3:
            atom3 = self.atoms[0][2]

        #Adicionar outras possibilidades de carregar cordenadas
        pkl_path   = os.path.join( _trajFolder, "")
        file_lists = glob.glob( pkl_path+"*.pkl" )
        self.bins  = len(file_lists)
        
        distance   = 0.0
        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )

        #...................................................................................
        if self.multipleDistance[0]:            
            with pymp.Parallel(self.nprocs) as p:
                for i in p.range( len(file_lists) ):
                    self.molecule.coordinates3 = ImportCoordinates3(file_lists[i])
                    dist_12 = self.molecule.coordinates3.Distance( atom1, atom2 )
                    dist_23 = self.molecule.coordinates3.Distance( atom2, atom3 )
                    distance = ( weight1 * dist_12) - ( weight2 * dist_23*-1)
                    rmodel = RestraintEnergyModel.Harmonic( distance, self.forceC )
                    restraint = RestraintMultipleDistance.WithOptions(energyModel = rmodel,  distances= [ [ atom2, atom1, weight1 ], [ atom2, atom3, weight2 ] ]) 
                
                    restraints['ReactionCoord'] = restraint

                    temp = file_lists[i][:-4]
                    temp = os.path.basename(temp)

                    mdfolder = os.path.join( self.baseName, temp )
                    mdRun = MD(self.molecule,mdfolder,self.mdMethod)
                    mdRun.ChangeDefaultParameters(self.mdParameters)
                    mdRun.RunProductionRestricted(self.equiNsteps,self.prodNsteps,_sample) 
        #......................................................................................
        else:
            with pymp.Parallel(self.nprocs) as p:
                for i in p.range( len(file_lists) ):
                    self.molecule.coordinates3 = ImportCoordinates3(file_lists[i])                    
                    distance    = self.molecule.coordinates3.Distance( atom1, atom2 )
                    rmodel      = RestraintEnergyModel.Harmonic( distance, self.forceC )
                    restraint   = RestraintDistance.WithOptions( energyModel = rmodel, point1=atom1, point2=atom2 ) 

                    restraints['ReactionCoord'] = restraint

                    temp = file_lists[i][:-4]
                    temp = os.path.basename(temp)
                    mdfolder = os.path.join( self.baseName, temp )
                    mdRun = MD(self.molecule,mdfolder,self.mdMethod)
                    mdRun.ChangeDefaultParameters(self.mdParameters)
                    mdRun.RunProductionRestricted(self.equiNsteps,self.prodNsteps,_sample) 

    #-----------------------------------------------------------------
    def RunTWODimensional(self,_trajFolder,_sample):
        '''
        Class method to execute two-dimesninal sampling 
        '''
        
        #Setting some local vars to ease the notation in the pDynamo methods
        #----------------------------------------------------------------
        atom1 = self.atoms[0][0]
        atom2 = self.atoms[0][1]
        atom3 = 0
        atom4 = 0
        atom5 = 0 
        atom6 = 0 

        weight1 = self.sigma_a3_a1[0]
        weight2 = self.sigma_a1_a3[0]
        weight3 = self.sigma_a3_a1[1]
        weight4 = self.sigma_a1_a3[1]

        if len(self.atoms[0]) == 3 and len(self.atoms[1]) == 2:
            atom3 = self.atoms[0][2]
            atom4 = self.atoms[1][0]
            atom5 = self.atoms[1][1]
        elif len(self.atoms[0]) == 3 and len(self.atoms[1]) == 3:
            atom3 = self.atoms[0][2]
            atom4 = self.atoms[1][0]
            atom5 = self.atoms[1][1]
            atom6 = self.atoms[1][2]
        elif len(self.atoms[0]) == 2 and len(self.atoms[1]) == 3:
            atom3 = self.atoms[1][0]
            atom4 = self.atoms[1][1]
            atom5 = self.atoms[1][2]
        elif len(self.atoms[0]) == 2 and len(self.atoms[1]) == 2:   
            atom3 = self.atoms[1][0]
            atom4 = self.atoms[1][1]
        else: print("Impossible combination!")

        #========================================================
        
        
        pkl_path   = os.path.join( _trajFolder, "")
        file_lists = glob.glob( pkl_path+"*.pkl" )
        self.bins  = len(file_lists)

        distance_1 = 0.0
        distance_2 = 0.0

        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )
        
        #............................................................
        if self.multipleDistance[0] and self.multipleDistance[1]:
            #Generate the initial structures for the row of the grid
            with pymp.Parallel(self.nprocs) as p:
                for i in p.range ( len(file_lists) ):  
                    
                    coordinate_file = file_lists[i]
                    temp = coordinate_file[:-4]
                    temp = os.path.basename(temp)
                    md_path = os.path.join(self.baseName, temp )
                    
                    if not os.path.exists( md_path ):              
                        self.molecule.coordinates3 = ImportCoordinates3( file_lists[i] )
                        distance_1  = self.molecule.coordinates3.Distance( atom1, atom2 ) - \
                                      self.molecule.coordinates3.Distance( atom2, atom3  )
                    
                        rmodel      =  RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                        restraint_1 =  RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = [ [ atom2, atom1, weight1 ],[ atom2, atom3, weight2 ] ] )
                        restraints["ReactionCoord"] = restraint_1
                                        
                        distance_2  = self.molecule.coordinates3.Distance( atom4, atom5 ) - \
                                      self.molecule.coordinates3.Distance( atom5, atom6 )

                        rmodel2     = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                        restraint_2 = RestraintMultipleDistance.WithOptions( energyModel = rmodel2,  distances = [ [ atom5, atom4, weight3 ],[ atom5, atom6, weight4 ] ] )
                        restraints["ReactionCoord2"] = restraint_2  

                        mdRun = MD(self.molecule,md_path,self.mdMethod)
                        mdRun.ChangeDefaultParameters(self.mdParameters)
                        mdRun.RunProductionRestricted(self.equiNsteps,self.prodNsteps,_sample)              
                    else:
                        continue
        #--------------------------------------------------------------------------------------------------------------
        elif self.multipleDistance[0] and not self.multipleDistance[1]:
            #Generate the initial structures for the row of the grid
            with pymp.Parallel(self.nprocs) as p:
                for i in p.range ( len(file_lists) ):  
                    coordinate_file = file_lists[i]
                    temp = coordinate_file[:-4]
                    temp = os.path.basename(temp)
                    md_path = os.path.join(self.baseName, temp)

                    if  not os.path.exists( md_path ):
                    #.----              
                        self.molecule.coordinates3 = ImportCoordinates3( file_lists[i] )
                        distance_1  = self.molecule.coordinates3.Distance( atom1, atom2 ) - \
                                      self.molecule.coordinates3.Distance( atom2, atom3 )
                    
                        rmodel      =  RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                        restraint_1 =  RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = [ [ atom2, atom1, weight1 ],[ atom2, atom3, weight2 ] ] )
                        restraints["ReactionCoord"] = restraint_1
                                        
                        distance_2  = self.molecule.coordinates3.Distance( atom4, atom5 )
                        rmodel2     = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                        restraint_2 = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom4, point2= atom5)
                                                                
                        restraints["ReactionCoord2"] = restraint_2  

                        mdRun = MD(self.molecule,md_path,self.mdMethod)
                        mdRun.ChangeDefaultParameters(self.mdParameters)
                        mdRun.RunProductionRestricted(self.equiNsteps,self.prodNsteps,_sample)              
                    
                    else:
                        continue
        #--------------------------------------------------------------------------------------------------------------
        elif not self.multipleDistance[0] and not self.multipleDistance[1]:
            with pymp.Parallel(self.nprocs) as p:
                for i in p.range ( len(file_lists) ):  
                    coordinate_file = file_lists[i]
                    temp = coordinate_file[:-4]
                    temp = os.path.basename(temp)
                    md_path = os.path.join(self.baseName,temp)
                    #.-----------------------------------------------------------------------------------
                    if  not os.path.exists( md_path ):
                                  
                        self.molecule.coordinates3 = ImportCoordinates3( file_lists[i] )
                       
                        distance_1  = self.molecule.coordinates3.Distance( atom1, atom2 )
                        rmodel      =  RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                        restraint_1 =  RestraintMultipleDistnce.WithOptions( energyModel = rmodel, point1= atom2, point2= atom1 )
                        restraints["ReactionCoord"] = restraint_1
                                        
                        distance_2  = self.molecule.coordinates3.Distance( self.atoms[1][0], self.atoms[1][1] )
                        rmodel2     = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                        restraint_2 = RestraintDistnce.WithOptions(energyModel = rmodel, point1= atom3, point2= atom4 )                      
                        restraints["ReactionCoord2"] = restraint_2  

                        mdRun = MD(self.molecule,self.mdParameters,mdfolder)

                        mdRun.RunProductionRestricted(self.equiSteps,self.prodSteps,_sample)              
                    
                    else:
                        continue

    #------------------------------------------------------------------------------------------------------------
    def Analysis(self):
        '''
        '''
        pass



#==================================================================================#
#================================END OF THE CLASS==================================#
#==================================================================================#           






