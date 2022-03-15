#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#FILE = UmbrellaSampling.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#==============================================================================

#-----------------------------------------------------
import os, glob, sys, shutil
#-----------------------------------------------------
from MolecularDynamics import MD 
#-----------------------------------------------------
import pymp
#-----------------------------------------------------
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
        self.multipleDistance   = [ False,False ]
        self.massConstraint     = True
        self.sigma_a1_a3        = [ 0.0,0.0 ]
        self.sigma_a3_a1        = [ 0.0,0.0 ]
        self.mdMethod           = mdMethod
        self.bins               = 0
        self.samplingFactor     = 0
        self.restart            = RESTART
        self.adaptative         = ADAPTATIVE
        self.angle              = False

        #qc/mm parameters 
        self.maxItQC            = 1000
        self.energyTolQC        = 1e-04
        self.densityTol         = 1e-08        
        self.timeStep           = 0.001
        self.mdParameters = { }        

    #====================================================================
    def ChangeDefaultParameters(self,_parameters):
        '''
        '''
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
                                'density_tolerancen':self.densityTol,\
                                'timeStep':self.timeStep,\
                                'energy_tolerance':self.energyTolQC
                            }

    #==========================================================================
    def SetMode(self,_RC):
        '''
        Class method to setup modes to be restrained
        '''
        #----------------------------------------------------------------------
        ndim = self.nDim # temp var to hold the index of the curren dim
        self.nDim += 1
        self.atoms.append(_RC.atoms)

        self.sigma_a1_a3[ndim]      = _RC.weight13
        self.sigma_a3_a1[ndim]      = _RC.weight31
        self.massConstraint         = _RC.massConstraint

        if len(_RC.atoms) == 3:
            self.multipleDistance[ndim] = True 
        if len(_RC.atoms) == 4:
            self.angle = True
    
    #============================================================================
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


    #============================================================================  
    def Run1DSampling(self,_trajFolder,_sample):
        '''
        Class method to execute one-dimensional sampling
        '''
        #-----------------------------------------------
        self.inputTraj      = _trajFolder
        self.samplingFactor = _sample
        #-----------------------------------------------               
        #Adicionar outras possibilidades de carregar cordenadas
        pkl_path        = os.path.join( _trajFolder, "")
        self.file_lists = glob.glob( pkl_path+"*.pkl" )
        self.bins       = len(self.file_lists)
        self.mdPaths    = []
        #------------------------------------------------------
        for i in range( len(self.file_lists) ):
            coordinate_file = self.file_lists[i]
            temp    = coordinate_file[:-4]
            temp    = os.path.basename(temp)
            md_path = os.path.join(self.baseName, temp )
            self.mdPaths.append(md_path)

        if self.restart:               
            for i in range( len(self.mdPaths) -1 , 0, -1 ):
                if os.path.exists( self.mdPaths[i] ):
                    self.mdPaths.remove( self.mdPaths[i] ) 
                    self.file_lists.remove( self.file_lists[i] ) 

        #-----------------------------------------------------
        if self.angle:
            self.Run1DSamplingDihedral()
        else:        
            if self.multipleDistance[0]:
                self.Run1DSamplingMultipleDistance()
            else:
                self.Run1DSamplingSimpleDistance()
    
    #==============================================================================
    def Run1DSamplingSimpleDistance(self):
        '''
        '''
        atom1 = self.atoms[0][0]
        atom2 = self.atoms[0][1]
        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )
        with pymp.Parallel(self.nprocs) as p:
            for i in p.range( len(self.file_lists) ):
                self.molecule.coordinates3 = ImportCoordinates3(self.file_lists[i])
                #------------------------------------------------------------                
                distance    = self.molecule.coordinates3.Distance( atom1, atom2 )
                rmodel      = RestraintEnergyModel.Harmonic( distance, self.forceC )
                restraint   = RestraintDistance.WithOptions( energyModel = rmodel, point1=atom1, point2=atom2 ) 
                restraints['M1'] = restraint
                #------------------------------------------------------------
                mdRun = MD(self.molecule,self.mdPaths[i],self.mdMethod)
                mdRun.ChangeDefaultParameters(self.mdParameters)
                mdRun.RunProductionRestricted(self.equiNsteps,self.prodNsteps,self.samplingFactor)

        self.molecule.DefineRestraintModel(None)

    #==============================================================================
    def Run1DSamplingMultipleDistance(self):
        '''
        '''
        atom1 = self.atoms[0][0]
        atom2 = self.atoms[0][1]
        atom3 = self.atoms[0][2]
        weight1 = self.sigma_a3_a1[0]
        weight2 = self.sigma_a1_a3[0] 
        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )
        with pymp.Parallel(self.nprocs) as p:
            for i in p.range( len(self.file_lists) ):
                self.molecule.coordinates3 = ImportCoordinates3(self.file_lists[i])
                #------------------------------------------------------------
                dist_12 = self.molecule.coordinates3.Distance( atom1, atom2 )
                dist_23 = self.molecule.coordinates3.Distance( atom2, atom3 )
                distance = ( weight1 * dist_12) - ( weight2 * dist_23*-1)
                #------------------------------------------------------------
                rmodel = RestraintEnergyModel.Harmonic( distance, self.forceC )
                restraint = RestraintMultipleDistance.WithOptions(energyModel = rmodel,  distances= [ [ atom2, atom1, weight1 ], [ atom2, atom3, weight2 ] ]) 
                restraints['M1'] = restraint
                #------------------------------------------------------------
                mdRun = MD(self.molecule,self.mdPaths[i],self.mdMethod)
                mdRun.ChangeDefaultParameters(self.mdParameters)
                mdRun.RunProductionRestricted(self.equiNsteps,self.prodNsteps,self.samplingFactor)
        
        self.molecule.DefineRestraintModel(None)
    #==============================================================================
    def Run1DSamplingDihedral(self):
        '''
        '''
        atom1 = self.atoms[0][0]
        atom2 = self.atoms[0][1]
        atom3 = self.atoms[0][2]
        atom4 = self.atoms[0][3]
        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )
        with pymp.Parallel(self.nprocs) as p:
            for i in p.range( len(self.file_lists) ):
                self.molecule.coordinates3 = ImportCoordinates3(self.file_lists[i])
                #------------------------------------------------------------
                angle = self.molecule.coordinates3.Dihedral( atom1, atom2, atom3, atom4 )
                #------------------------------------------------------------
                rmodel      = RestraintEnergyModel.Harmonic( angle, self.forceC )
                restraint   = RestraintMultipleDistance.WithOptions(energyModel = rmodel, point1=atom1,
                                                                                          point2=atom2,
                                                                                          point3=atom3,
                                                                                          point4=atom4) 
                restraints['M1'] = restraint
                #------------------------------------------------------------
                mdRun = MD(self.molecule,self.mdPaths[i],self.mdMethod)
                mdRun.ChangeDefaultParameters(self.mdParameters)
                mdRun.RunProductionRestricted(self.equiNsteps,self.prodNsteps,self.samplingFactor)
        #---------------------------------------
        self.molecule.DefineRestraintModel(None)

    #==============================================================================
    def Run2DSampling(self,_trajFolder,_sample):
        '''
        Class method to execute two-dimesninal sampling 
        ''' 
        self.samplingFactor = _sample
        pkl_path            = os.path.join( _trajFolder, "")
        self.file_lists     = glob.glob( pkl_path+"*.pkl" )
        self.mdPaths        = []
        #-----------------------------------------------
        for i in range( len(self.file_lists) ):
            coordinate_file = self.file_lists[i]
            temp    = coordinate_file[:-4]
            temp    = os.path.basename(temp)
            md_path = os.path.join(self.baseName, temp )
            self.mdPaths.append(md_path)
        
        if self.restart:               
            for i in range(len(self.mdPaths)-1,0,-1 ):
                if os.path.exists( self.mdPaths[i] ):
                    self.mdPaths.remove( self.mdPaths[i] ) 
                    self.file_lists.remove( self.file_lists[i] )   

        self.bins = len(self.file_lists)

        #-----------------------------------------------
        self.EnergyRef = self.molecule.Energy()
        self.forceCRef = self.forceC
        #-----------------------------------------------
        if self.angle: 
            self.Run2DSamplingDihedral()
        else:
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
        #-----------------------
        atom1 = self.atoms[0][0]
        atom2 = self.atoms[0][1]
        atom3 = self.atoms[0][2]
        atom4 = self.atoms[1][0]
        atom5 = self.atoms[1][1]
        atom6 = self.atoms[1][2]
        #-----------------------
        weight1 = self.sigma_a3_a1[0]
        weight2 = self.sigma_a1_a3[0]
        weight3 = self.sigma_a3_a1[1]
        weight4 = self.sigma_a1_a3[1]
        #-----------------------
        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )
        #-------------------------------------
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
                restraint   =  RestraintMultipleDistance.WithOptions(energyModel = rmodel, distances = [ [ atom2, atom1, weight1 ],[ atom2, atom3, weight2 ] ] )
                restraints["M1"] = restraint
                #--------------------------------------------------------               
                dist45      = self.molecule.coordinates3.Distance(atom4, atom5)
                dist56      = self.molecule.coordinates3.Distance(atom5, atom6)
                distance_2  = ( weight1 * dist45) - ( weight2 * dist56*-1)
                #--------------------------------------------------------
                rmodel2     = RestraintEnergyModel.Harmonic(distance_2,self.forceC)
                restraint   = RestraintMultipleDistance.WithOptions(energyModel = rmodel2,distances = [ [ atom5, atom4, weight3 ],[ atom5, atom6, weight4 ] ] )
                restraints["M2"] = restraint 
                #--------------------------------------------------------
                if self.adaptative:
                    self.ChangeConvergenceParameters()
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
                if self.adaptative:
                    self.ChangeDefaultParameters()
                mdRun = MD(self.molecule,self.mdPaths[i],self.mdMethod)
                mdRun.ChangeDefaultParameters(self.mdParameters)
                mdRun.RunProductionRestricted(self.equiNsteps,self.prodNsteps,self.samplingFactor)
        #---------------------------------------
        self.molecule.DefineRestraintModel(None)
        #.......................................

    #==========================================================================================
    def Run2DSimpleDistance(self):
        '''
        Run sampling with the two dimensions set as simple distances types.
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
                distance_1       = self.molecule.coordinates3.Distance(atom1, atom2)                
                rmodel           = RestraintEnergyModel.Harmonic(distance_1, self.forceC)
                restraint        = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom1, point2= atom2)
                restraints["M1"] = restraint
                #-----------------------------------------------------------------------         
                distance_2       = self.molecule.coordinates3.Distance(atom3, atom4)
                rmodel           = RestraintEnergyModel.Harmonic(distance_2, self.forceC)
                restraint        = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom3, point2= atom4)
                restraints["M2"] = restraint  
                #-----------------------------------------------------------------------
                if self.adaptative:
                    self.ChangeConvergenceParameters()
                #-----------------------------------------------------------------------
                mdRun = MD(self.molecule,self.mdPaths[i],self.mdMethod)
                mdRun.ChangeDefaultParameters(self.mdParameters)
                mdRun.RunProductionRestricted(self.equiNsteps,self.prodNsteps,self.samplingFactor)
        #---------------------------------------
        self.molecule.DefineRestraintModel(None)
        #.......................................
        
    #===========================================================================================
    def Run2DSamplingDihedral(self):
        '''
        '''
        atom1 = self.atoms[0][0]
        atom2 = self.atoms[0][1]
        atom3 = self.atoms[0][2]
        atom4 = self.atoms[0][3]
        atom5 = self.atoms[1][0]
        atom6 = self.atoms[1][1]
        atom7 = self.atoms[1][2]
        atom8 = self.atoms[1][3]

        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )
        with pymp.Parallel(self.nprocs) as p:
            for i in p.range( len(self.file_lists) ):
                self.molecule.coordinates3 = ImportCoordinates3(self.file_lists[i])
                #------------------------------------------------------------
                angle_1 = self.molecule.coordinates3.Dihedral( atom1, atom2, atom3, atom4 )
                #------------------------------------------------------------
                rmodel      = RestraintEnergyModel.Harmonic( angle_1, self.forceC )
                restraint   = RestraintMultipleDistance.WithOptions(energyModel = rmodel, point1=atom1,
                                                                                          point2=atom2,
                                                                                          point3=atom3,
                                                                                          point4=atom4) 
                restraints['M1'] = restraint
                angle_2 = self.molecule.coordinates3.Dihedral( atom5, atom6, atom7, atom8 )
                #------------------------------------------------------------
                rmodel      = RestraintEnergyModel.Harmonic( angle_2, self.forceC )
                restraint   = RestraintMultipleDistance.WithOptions(energyModel = rmodel, point1=atom5,
                                                                                          point2=atom6,
                                                                                          point3=atom7,
                                                                                          point4=atom8) 
                restraints['M2'] = restraint
                #------------------------------------------------------------
                mdRun = MD(self.molecule,self.mdPaths[i],self.mdMethod)
                mdRun.ChangeDefaultParameters(self.mdParameters)
                mdRun.RunProductionRestricted(self.equiNsteps,self.prodNsteps,self.samplingFactor)
        #---------------------------------------
        self.molecule.DefineRestraintModel(None)

    #===========================================================================================
    def Finalize(self):
        '''
        Reorganize frames and concatenate in a single trajectory folder 
        '''
        fsize = int(self.prodNsteps/self.samplingFactor)
        
        self.concFolder = os.path.join(self.baseName,"concatenated_trajectory.ptGeo")
        
        if not os.path.exists(self.concFolder):
            os.makedirs( self.concFolder )

        if self.nDim == 1:
            pkl_path = self.baseName + "/frame*/production*/frame*.pkl"
            pkl_paths = glob.glob( pkl_path )
            pkl_paths.sort() 
            
            newNames = []
            cnt = 0
            for i in range(self.bins):
                for j in range(fsize):
                    newName = os.path.join(self.concFolder,"frame{}.pkl".format( j + i*fsize ) )
                    shutil.copy(pkl_paths[cnt],newName)
                    cnt+=1

            Duplicate(self.concFolder,self.baseName+".dcd",self.molecule)

#==================================================================================#
#================================END OF THE CLASS==================================#
#==================================================================================#           






