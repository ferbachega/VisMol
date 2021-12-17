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
        self.nprocs             = 1
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
        self.samplingFactor     = self.prodnsteps/20 # sampling 20 frames per window

        #qc/mm parameters 
        self.maxItQC            = 1000
        self.energyTolQC        = 1e-04
        self.densityTol         = 1e-08        
        self.timeStep           = 0.001
        self.mdParameters = { }        


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
        self.sigma_a1_a3[ndim] = 1 
        self.sigma_a3_a1[ndim] = -1

        if len(_atoms) == 3:
            self.multipleDistance = True

        #.---------------------------
        if self.multipleDistance[0]:
            #.----------------------
            if self.massConstraint:
                atomic_n1 = self.molecule.atoms.items[ self.atoms[ndim][0] ].atomicNumber
                atomic_n3 = self.molecule.atoms.items[ self.atoms[ndim][2] ].atomicNumber
                mass_a1 = atomic_mass[ atomic_n1 ]
                mass_a3 = atomic_mass[ atomic_n3 ]
                self.sigma_a1_a3[ndim] = mass_a1 /(mass_a1+mass_a3)
                self.sigma_a3_a1[ndim] = mass_a3 /(mass_a1+mass_a3)
                self.sigma_a3_a1[ndim] = self.sigma_a3_a1[ndim]*-1

                
    #----------------------------------------------------------------------------------  
    def RunONEDimensional(self,_trajFolder,_sample):
        '''
        Class method to execute one-dimensional sampling
        '''
        if not _sample:
            self.samplingFactor = 0
        #Adicionar outras possibilidades de carregar cordenadas
        pkl_path   = os.path.join( _trajFolder, "")
        file_lists = glob.glob( pkl_path+"*.pkl" )

        distance = 0.0
        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )

        self.molecule.energyModel.qcModel.converger.SetOptions( maximumSCFCycles = self.maxItQC )

        #...................................................................................
        if self.multipleDistance[0]:
            for i in range( len(file_lists) ):
                self.molecule.coordinates3 = Unpickle(file_lists[i])
                distance = self.molecule.coordinates3( self.atoms[0][0], self.atoms[0][1] )
                rmodel = RestraintEnergyModel.Harmonic( distance, self.forceC )
                restraint = RestraintMultipleDistnce.WithOptions(energyModel = rmodel, [\
                                                                [ point1= self.atoms[0][1], point1= self.atoms[0][0], weight=self.sigma_a1_a3[0] ],\
                                                                [ point1= self.atoms[0][1], point1= self.atoms[0][2], weight=self.sigma_a3_a1[0] ]\
                                                                ] ) 

                restraints['ReactionCoord'] = restraint
                mdfolder = os.path.join( self.baseName, file_lists[i]+"_MD" )
                mdRun = MD(self.molecule,mdfolder,self.mdMethod)
                MD.ChengeDefaultParameters(self.mdParameters)
                MD.RunProductionRestricted(self.equiSteps,self.prodSteps,self.mdParameters['sampling_factor']) 
        #......................................................................................
        else:
            for i in range( len(file_lists) ):
                self.molecule.coordinates3 = Unpickle(file_lists[i])
                distance = self.molecule.coordinates3( self.atoms[0][0], self.atoms[0][1] )
                rmodel = RestraintEnergyModel.Harmonic( distance, self.forceC )
                restraint = RestraintDistnce.WithOptions(energyModel = rmodel, point1= self.atoms[0][0], point2= self.atoms[0][1] ) 

                restraints['ReactionCoord'] = restraint
                mdfolder = os.path.join( self.baseName, file_lists[i]+"_MD")
                mdRun = MD(self.molecule,mdfolder,self.mdMethod)
                mdRun.ChengeDefaultParameters(self.mdParameters)
                MD.RunProductionRestricted(self.equiNsteps,self.prodNsteps,self.mdParameters['sampling_factor']) 

    #-----------------------------------------------------------------
    def RunTWODimensional(self,_trajFolder):
        '''
        Class method to execute two-dimesninal sampling 
        '''
        pkl_path   = os.path.join( _trajFolder, "")
        file_lists = glob.glob( pkl_path+"*.pkl" )

        distance_1 = 0.0
        distance_2 = 0.0

        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )

        self.molecule.energyModel.qcModel.converger.SetOptions( maximumSCFCycles = self.maxItQC )
        
        #............................................................
        if self.multipleDistance[0] and self.multipleDistance[1]:
            #Generate the initial structures for the row of the grid
            with pymp.Parallel(self.nprocs) as p:
                for i in p.range ( len(file_lists) ):  
                    coordinate_file = files_list[i]
                    temp = coordinate_file[:-4]
                    temp = os.path.basename(temp)
                    md_path = os.path.join(self.baseName, file_lists[i]+"_MD")

                    if  not os.path.exists( md_path ):
                    #.----              
                        self.molecule.coordinates3 = Unpickle( file_lists[i] )
                        distance_1  = self.molecule.coordinates3.Distance( self.atoms[0][0], self.atoms[0][1] ) - \
                                      self.molecule.coordinates3.Distance( self.atoms[0][1], self.atoms[0][2] )
                    
                        rmodel      =  RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                        restraint_1 =  RestraintMultipleDistnce.WithOptions( energyModel = rmodel, [\
                                                                [ point1= self.atoms[0][1], point1= self.atoms[0][0], weight=self.sigma_a1_a3[0] ],\
                                                                [ point1= self.atoms[0][1], point1= self.atoms[0][2], weight=self.sigma_a3_a1[0] ] \
                                                                ] )
                        restraints["ReactionCoord"] = restraint_1
                    
                    
                        distance_2  = self.molecule.coordinates3.Distance( self.atoms[1][0], self.atoms[1][1] ) - \
                                      self.molecule.coordinates3.Distance( self.atoms[1][1], self.atoms[1][2] )

                        rmodel2     = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                        restraint_1 = RestraintMultipleDistnce.WithOptions(energyModel = rmodel, [\
                                                                [ point1= self.atoms[1][1], point1= self.atoms[1][0], weight=self.sigma_a1_a3[1] ],\
                                                                [ point1= self.atoms[1][1], point1= self.atoms[1][2], weight=self.sigma_a3_a1[1] ] \
                                                                ] )
                        
                        restraints["ReactionCoord2"] = restraint_2  
                        mdRun = MD(self.molecule,mdfolder,self.mdMethod)
                        mdRun.ChengeDefaultParameters(self.mdParameters)
                        mdRun.RunProductionRestricted(self.equiSteps,self.prodSteps,self.mdParameters['sampling_factor'])              
                    else:
                        continue
        #--------------------------------------------------------------------------------------------------------------
        elif self.multipleDistance[0] and not self.multipleDistance[1]:
            #Generate the initial structures for the row of the grid
            with pymp.Parallel(self.nprocs) as p:
                for i in p.range ( len(file_lists) ):  
                    coordinate_file = files_list[i]
                    temp = coordinate_file[:-4]
                    temp = os.path.basename(temp)
                    md_path = os.path.join(self.baseName, file_lists[i]+"_MD")

                    if  not os.path.exists( md_path ):
                    #.----              
                        self.molecule.coordinates3 = Unpickle( file_lists[i] )
                        distance_1  = self.molecule.coordinates3.Distance( self.atoms[0][0], self.atoms[0][1] ) - \
                                      self.molecule.coordinates3.Distance( self.atoms[0][1], self.atoms[0][2] )
                    
                        rmodel      =  RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                        restraint_1 =  RestraintMultipleDistnce.WithOptions( energyModel = rmodel, [\
                                                                [ point1= self.atoms[0][1], point1= self.atoms[0][0], weight=self.sigma_a1_a3[0] ],\
                                                                [ point1= self.atoms[0][1], point1= self.atoms[0][2], weight=self.sigma_a3_a1[0] ] \
                                                                ] )
                        restraints["ReactionCoord"] = restraint_1
                    
                    
                        distance_2  = self.molecule.coordinates3.Distance( self.atoms[1][0], self.atoms[1][1] )
                        rmodel2     = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                        restraint_1 = RestraintDistnce.WithOptions(energyModel = rmodel, point1= self.atoms[1][1], point1= self.atoms[1][0] )
                                                                
                        restraints["ReactionCoord2"] = restraint_2  
                        mdRun = MD(self.molecule,mdfolder,self.mdMethod)
                        mdRun.ChengeDefaultParameters(self.mdParameters)
                        mdRun.RunProductionRestricted(self.equiSteps,self.prodSteps,self.mdParameters['sampling_factor'])              
                    
                    else:
                        continue
        #--------------------------------------------------------------------------------------------------------------
        elif not self.multipleDistance[0] and not self.multipleDistance[1]:
            with pymp.Parallel(self.nprocs) as p:
                for i in p.range ( len(file_lists) ):  
                    coordinate_file = files_list[i]
                    temp = coordinate_file[:-4]
                    temp = os.path.basename(temp)
                    md_path = os.path.join(self.baseName, file_lists[i]+"_MD")
                    #.-----------------------------------------------------------------------------------
                    if  not os.path.exists( md_path ):
                                  
                        self.molecule.coordinates3 = Unpickle( file_lists[i] )
                        distance_1  = self.molecule.coordinates3.Distance( self.atoms[0][0], self.atoms[0][1] )
                    
                        rmodel      =  RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                        restraint_1 =  RestraintMultipleDistnce.WithOptions( energyModel = rmodel, point1= self.atoms[0][1], point1= self.atoms[0][0] )
                
                        restraints["ReactionCoord"] = restraint_1
                    
                    
                        distance_2  = self.molecule.coordinates3.Distance( self.atoms[1][0], self.atoms[1][1] )
                        rmodel2     = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                        restraint_1 = RestraintDistnce.WithOptions(energyModel = rmodel, point1= self.atoms[1][1], point1= self.atoms[1][0] )
                                                                
                        restraints["ReactionCoord2"] = restraint_2  
                        mdRun = MD(self.molecule,self.mdParameters,mdfolder)
                        MD.RunProductionRestricted(self.equiSteps,self.prodSteps,self.mdParameters['sampling_factor'])              
                    
                    else:
                        continue



#==================================================================================#
#================================END OF THE CLASS==================================#
#==================================================================================#           






