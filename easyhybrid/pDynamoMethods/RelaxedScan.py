#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#FILE = RelaxedScan.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#==============================================================================

import pymp
import GeometrySearcher 

from commonFunctions import *
from pMolecule import *
from GeometrySearcher import* 

#*****************************************************************************
class SCAN:
    '''
    Class to setup and execute relaxed surface scan procedure
    '''
    
    def __init__(self,_system,_baseFolder,_optimizer):
        '''
        Class constructor
        '''
        self.baseName           = _baseFolder
        self.molecule           = _system 
        self.nDim               = 0
        self.reactionCoordinate1= []
        self.reactionCoordinate2= []
        self.atoms              = [] 
        self.nprocs             = 1
        self.textLog            = " "
        self.energies           = []
        self.DMINIMUM           = [ 0.0, 0.0 ]
        self.DINCREMENT         = [ 0.0, 0.0 ]
        self.forceC             = 2500.0
        self.massConstraint     = True
        self.multipleDistance   = [False,False]
        self.nsteps             = [ 1, 1 ]
        self.logFreq            = 50
        self.maxIt              = 300
        self.maxItQC            = 200
        self.rmsGT              = 0.1
        self.optmizer           = _optimizer
        self.sigma_a1_a3        = [0.0,0.0]
        self.sigma_a3_a1        = [0.0,0.0]
        self.real_distance_1    = []
        self.real_distance_2    = []
        self.text = ""
    
        self.baseName = os.path.join( self.baseName, "Relaxed_Surface_Scan")
        if not os.path.exists( self.baseName ):
            os.makedirs( self.baseName )
        if not os.path.exists(  os.path.join( self.baseName, "ScanTraj.ptGeo" ) ):
            os.makedirs(  os.path.join( self.baseName, "ScanTraj.ptGeo" ) )

        #set the parameters dict for the geometry search classes


        self.GeoOptPars =   { 
                                "method": self.optmizer           ,\
                                "logFrequency": self.logFreq      ,\
                                "savePdb": "false"                ,\
                                "saveTraj": "false"               ,\
                                "maxIterations":self.maxIt        ,\
                                "maxIrerations_QC": self.maxItQC  ,\
                                "rmsGradient": self.rmsGT
                            }
    #-------------------------------------------------------------------------------------
    def ChangeDefaultParameters(self,_parameters):
        '''
        Class method to alter deafult parameters
        '''
        if "rmsGradient" in _parameters:
            self.rmsGT = _parameters['rmsGradient']
        if "maxIrerations" in _parameters:
            self.maxIt = _parameters['maxIterations']
        if "maxIterations_QC" in _parameters:
            self.maxItQC
        if "log_frequency" in _parameters:
            self.logFreq = _paremeters["log_frequency"]

        self.GeoOptPars =   { 
                            "method": self.optmizer           ,\
                            "logFrequency": self.logFreq      ,\
                            "not_save_pdb": "false"           ,\
                            "saveTraj": "false"               ,\
                            "maxIterations":self.maxIt        ,\
                            "maxIrerations_QC": self.maxItQC  ,\
                            "rmsGradient": self.rmsGT
                            }
   

    #---------------------------------------------------------------------------------------
    def SetReactionCoord(self,_atoms,_dincre,_massConstraint):
        '''
        Set reaction coordinate, determining initial parameters from the atoms information
        '''
        self.massConstraint = _massConstraint
        ndim = self.nDim # temp var to hold the index of the curren dim
        self.nDim += 1
        self.atoms.append(_atoms)
        self.DINCREMENT[ndim] = _dincre
        self.sigma_a1_a3[ndim] = 1.0
        self.sigma_a3_a1[ndim] = -1.0

        if len(_atoms) == 3:
            self.multipleDistance[ndim] = True

        #.---------------------------
        if self.multipleDistance[0]:
            #.----------------------
            if self.massConstraint:
                atomic_n1 = self.molecule.atoms.items[ self.atoms[ndim][0] ].atomicNumber
                atomic_n3 = self.molecule.atoms.items[ self.atoms[ndim][2] ].atomicNumber
                mass_a1 = GetAtomicMass(atomic_n1)
                mass_a3 = GetAtomicMass(atomic_n3)
                self.sigma_a1_a3[ndim] = mass_a1 /(mass_a1+mass_a3)
                self.sigma_a3_a1[ndim] = mass_a3 /(mass_a1+mass_a3)
                self.sigma_a3_a1[ndim] = self.sigma_a3_a1[ndim]*-1
                dist_a1_a2 = self.molecule.coordinates3.Distance( self.atoms[ndim][0], self.atoms[ndim][1] )
                dist_a2_a3 = self.molecule.coordinates3.Distance( self.atoms[ndim][1], self.atoms[ndim][2] )
                self.DMINIMUM[ndim] = ( self.sigma_a1_a3[ndim] * dist_a1_a2) - ( self.sigma_a3_a1[ndim] * dist_a2_a3*-1)

            #.----------------------
            else:
                dist_a1_a2 = self.molecule.coordinates3.Distance( self.atoms[ndim][0], self.atoms[ndim][1] )
                dist_a2_a3 = self.molecule.coordinates3.Distance( self.atoms[ndim][1], self.atoms[ndim][2] )
                self.DMINIMUM[ndim] =  dist_a1_a2 - dist_a2_a3
        #.--------------------------       
        else:
            self.DMINIMUM[ndim] = self.molecule.coordinates3.Distance( self.atoms[ndim][0], self.atoms[ndim][1] ) 


    #--------------------------------------------------------------------------------------------------
    def RunONEDimensionSCAN(self,_nsteps):
        '''
        Execute the relaxed scan with one reaction coordinate
        '''

        #Setting some local vars to ease the notation in the pDynamo methods
        #----------------------------------
        atom1 = self.atoms[0][0]
        atom2 = self.atoms[0][1]
        atom3 = 0

        weight1 = self.sigma_a3_a1[0]
        weight2 = self.sigma_a1_a3[0]

        if len(self.atoms[0]) == 3:
            atom3 = self.atoms[0][2]
        #---------------------------------

        self.text += "x RC1 Energy\n" 

        restraints = RestraintModel ( )
        self.molecule.DefineRestraintModel( restraints )
        distance = 0.0

        en0 = self.molecule.Energy()

        #---------------------------------
        if self.multipleDistance[0]:
            for i in range(_nsteps):     
                distance = self.DMINIMUM[0] + self.DINCREMENT[0]* float( i )
                rmodel = RestraintEnergyModel.Harmonic( distance, self.forceC )
                restraint = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom2, atom1, weight1 ], [ atom2, atom3, weight2 ] ] )
                restraints['ReactionCoord'] =  restraint            

                relaxRun = GeometrySearcher(self.molecule, self.baseName  )
                relaxRun.ChangeDefaultParameters(self.GeoOptPars)
                relaxRun.Minimization(self.optmizer)
                
                self.energies.append( self.molecule.Energy() )
                self.reactionCoordinate1.append( self.molecule.coordinates3.Distance( atom1 , atom2  ) - self.molecule.coordinates3.Distance( atom2, atom3  ) ) 
                Pickle( os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}.pkl".format(i) ), self.molecule.coordinates3 )
        else:
             for i in range(_nsteps):       
                distance = self.DMINIMUM[0] + self.DINCREMENT[0] * float( i )
                rmodel = RestraintEnergyModel.Harmonic( distance, self.forceC )
                restraint = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom1, point2= atom2 )
                restraints['ReactionCoord'] =  restraint            

                relaxRun = GeometrySearcher(self.molecule,"test")
                relaxRun.ChangeDefaultParameters(self.GeoOptPars)
                relaxRun.Minimization(self.optmizer)

                self.energies.append( self.molecule.Energy() )
                self.reactionCoordinate1.append( self.molecule.coordinates3.Distance( atom1 , atom2  ) )  
                Pickle( os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}.pkl".format(i) ), self.molecule.coordinates3 )   
        
        print(self.reactionCoordinate1[0])
        input()

        for i in range(_nsteps):
            self.text += "{} {} {} \n".format( i,self.reactionCoordinate1[i],self.energies[i]-en0 )
        
        trajNameDCD = self.baseName + ".dcd"
        trajName = os.path.join( self.baseName, "ScanTraj.ptGeo" )
        Duplicate( trajName, trajNameDCD, self.molecule )
        #clean the restraints assigned to the sytem
        textLog = open( self.baseName+"_scan1D.log", "w" ) 
        textLog.write( self.text )
        textLog.close()
        self.molecule.DefineRestraintModel ( None )

    #---------------------------------------------------------------------------------------------------
    def RunTWODimensionSCAN(self,_nsteps_x,_nsteps_y):
        '''
        Execute the relaxed scan with two reaction coordinate
        '''
        #Setting some local vars to ease the notation in the pDynamo methods
        #----------------------------------
        atom1 = self.atoms[0]
        atom2 = self.atoms[1]
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

        #---------------------------------

        self.text += "x y RC1 RC2 Energy\n" 

        restraints = RestraintModel ( )
        self.molecule.DefineRestraintModel( restraints )
        rmodel = RestraintEnergyModel.Harmonic( )
        distance_1 = 0.0 
        distance_2 = 0.0

        N = _nsteps_x
        M = _nsteps_y

        
        #Define the origin point of the relaxed surface scan, AKA the 0,0 point
        coordinate_file = os.path.join( self.baseName ,"ScanTraj","frame{}_{}.pkl".format( 0, 0 ) )
        relaxRun0 = Geometry_Optimization( self.molecule, self.GeoOptPars, self.baseName )
        relaxRun0.Minimization()
        en0 = self.molecule.Energy()
        Pickle( coordinate_file, self.molecule.coordinates3 )
        
        #............................................................
        if self.multipleDistance[0] and self.multipleDistance[1]:
            #Generate the initial structures for the row of the grid
            with pymp.Parallel(self.nprocs) as p:
                for i in p.range ( 1, M ):  
                    #.----              
                    distance_1 = self.DINCREMENT[0] * float(i) + self.DMINIMUM[0]
                    rmodel  =  RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                    restraint_1=  RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = [ [ atom2, atom1, weight1 ] , [ atom2, atom3, weight2 ] ] )
                    restraints["ReactionCoord"] = restraint_1
                    
                    j=0

                    distance_2  = self.DINCREMENT[1] * float(j) + self.DMINIMUM[1]
                    rmodel2      = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                    restraint_2 = RestraintMultipleDistance.WithOptions( energyModel = rmodel2, distances = [ [ atom5, atom4, weight3 ],[ atom5, atom6, weight4 ] ] )

                    restraints["ReactionCoord2"] = restraint_2  
                                        
                    initCoordinateFile = os.path.join( self.baseName,"ScanTraj", "frame{}_{}.pkl".format(0,0) ) 
                    self.molecule.coordinates3 = Unpickle( initCoordinateFile )             
                    coordinateFile = os.path.join( self.baseName,"ScanTraj", "frame{}_{}.pkl".format( i, j ) )
                    relaxRun = Geometry_Optimization( self.molecule, self.baseName )
                    relaxRun.ChangeDefaultParameters(self.GeoOptPars)
                    relaxRun.Minimization()
                    Pickle( coordinateFile, self.molecule.coordinates3 ) 
            
            #...........................................................................
            with pymp.Parallel(self.nprocs) as p:
                #Pergomr the calculations for the rest of the grid
                for i in p.range ( 0, M ):
                    distance_1 = self.DINCREMENT[0] * float(i) + self.DMINIMUM[0]
                    rmodel  =  RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                    restraint_1=  RestraintMultipleDistance.WithOptions(energyModel = rmodel, distances= [ [ atom2, atom1, weight1 ],[ atom2, atom3, weight2 ] ] )
                    restraints["ReactionCoord"] = restraint_1

                    for j in range( 1,N ):
                        distance_2  = self.DINCREMENT[1] * float(j) + self.DMINIMUM[1]
                        rmodel2     = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                        restraint_2 = RestraintMultipleDistance.WithOptions(energyModel = rmodel2, distances = [ [ atom5, atom4, weight3 ],[ atom5, atom6, weight4 ] ] )
                        restraints["ReactionCoord2"] = restraint_2  
                        
                        initCoordinateFile = ""
                        if  j==0:
                            initCoordinateFile = os.path.join( self.baseName,"ScanTraj" , "frame{}_{}.pkl".format( i, 0 ) )
                        elif j>0:
                            initCoordinateFile = os.path.join( self.baseName,"ScanTraj" , "frame{}_{}.pkl".format( i, j-1 ) )                 
                        
                        self.molecule.coordinates3 = Unpickle( initCoordinateFile )             
                        coordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( i, j ) )
                        relaxRun = Geometry_Optimization( self.molecule, self.baseName  )
                        relaxRun.ChangeDefaultParameters(self.GeoOptPars)
                        relaxRun.Minimization()
                        Pickle( coordinateFile, self.molecule.coordinates3 )


            #Recuperate energy and reaction coordinate information   
            for i in range(M):
                for j in range(N):
                    molecule = self.molecule
                    molecule.coordinates3 = Unpickle( os.path.join( self.folder_trj , "frame{}_{}.pkl".format( i, j ) ) )
                    Energy = molecule.Energy()
                    RD1 =   molecule.coordinates3.Distance( atom1, atom2 ) - molecule.coordinates3.Distance( atom2, atom3 )   
                    RD2 =   molecule.coordinates3.Distance( atom4, atom5 ) - molecule.coordinates3.Distance( atom5, atom6 )
                    self.text += "{} {} {} {} {}\n".format( i,j,RD1,RD2,Energy-en0 )
        
        #---------------------------------------------------------------------------------------------------------------------------
        elif self.multipleDistance[0] and self.multipleDistance[1] == False:
            with pymp.Parallel(self.nprocs) as p:
                for i in p.range ( 1, M ):  
                    #.----              
                    distance_1 = self.DINCREMENT[0] * float(i) + self.DMINIMUM[0]
                    rmodel  =  RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                    restraint_1=  RestraintMultipleDistance.WithOptions(energyModel = rmodel, distances = [ [ atom2, atom1, weight1 ],[ atom2, atom3, weight2 ] ] )
                    restraints["ReactionCoord"] = restraint_1
                    
                    j=0

                    distance_2  = self.DINCREMENT[1] * float(j) + self.DMINIMUM[1]
                    rmodel2      = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                    restraint_2 = RestraintDistance.WithOptions(energyModel = rmodel2, point1=atom4, point2=atom5 )

                    restraints["ReactionCoord2"] = restraint_2  
                                        
                    initCoordinateFile = os.path.join( self.baseName,"scan_trj", "frame{}_{}.pkl".format(0,0) ) 
                    self.molecule.coordinates3 = Unpickle( initCoordinateFile )             
                    coordinateFile = os.path.join( self.baseName,"scan_trj", "frame{}_{}.pkl".format( i, j ) )
                    relaxRun = Geometry_Optimization( self.molecule, self.baseName )
                    relaxRun.Minimization()
                    Pickle( coordinateFile, self.molecule.coordinates3 ) 
            
            #...........................................................................
            with pymp.Parallel(self.nprocs) as p:
                #Pergomr the calculations for the rest of the grid
                for i in p.range ( 0, M ):
                    distance_1 = self.DINCREMENT[0] * float(i) + self.DMINIMUM[0]
                    rmodel  =  RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                    restraint_1=  RestraintMultipleDistance.WithOptions(energyModel = rmodel, distances= [ [ atom2, atom1, weight1 ],[ atom2, atom3, weight2 ] ] )
                    restraints["ReactionCoord"] = restraint_1

                    for j in range( 1,N ):
                        distance_2  = self.DINCREMENT[1] * float(j) + self.DMINIMUM[1]
                        rmodel2     = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                        restraint_2 = RestraintMultipleDistance.WithOptions(energyModel = rmodel2, distances = [ [ atom5, atom4, weight3 ],[ atom5, atom6, weight4 ] ] )
                        restraints["ReactionCoord2"] = restraint_2  
                        
                        initCoordinateFile = ""
                        if  j==0:
                            initCoordinateFile = os.path.join( self.baseName,"scan_trj" , "frame{}_{}.pkl".format( i, 0 ) )
                        elif j>0:
                            initCoordinateFile = os.path.join( self.baseName,"scan_trj" , "frame{}_{}.pkl".format( i, j-1 ) )                 
                        
                        self.molecule.coordinates3 = Unpickle( initCoordinateFile )             
                        coordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( i, j ) )
                        relaxRun = Geometry_Optimization( self.molecule, self.baseName  )
                        relaxRun.ChangeDefaultParameters( self.GeoOptPars )
                        relaxRun.Minimization()
                        Pickle( coordinateFile, self.molecule.coordinates3 )


            #Recuperate energy and reaction coordinate information   
            for i in range(M):
                for j in range(N):
                    molecule = self.molecule
                    molecule.coordinates3 = Unpickle( os.path.join( self.folder_trj , "frame{}_{}.pkl".format( i, j ) ) )
                    Energy = molecule.Energy()
                    RD1 =   molecule.coordinates3.Distance( atom1, atom2 ) - molecule.coordinates3.Distance( atom2, atom3 )   
                    RD2 =   molecule.coordinates3.Distance( atom4, atom5 ) - molecule.coordinates3.Distance( atom5, atom6 )
                    self.text += "{} {} {} {} {}\n".format( i,j,RD1,RD2,Energy-en0 )


        elif self.multipleDistance[0] == False and self.multipleDistance[1] == False:
            with pymp.Parallel(self.nprocs) as p:
                for i in p.range ( 1, M ):  
                    #.----              
                    distance_1 = self.DINCREMENT[0] * float(i) + self.DMINIMUM[0]
                    rmodel  =  RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                    restraint_1=  RestraintDistance.WithOptions(energyModel = rmodel,  point1=atom1, point2=atom2  )
                    restraints["ReactionCoord"] = restraint_1
                    
                    j=0

                    distance_2  = self.DINCREMENT[1] * float(j) + self.DMINIMUM[1]
                    rmodel2      = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                    restraint_2 = RestraintDistance.WithOptions(energyModel = rmodel2, point1=atom3, point2=atom4 )

                    restraints["ReactionCoord2"] = restraint_2  
                                        
                    initCoordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format(0,0) ) 
                    self.molecule.coordinates3 = Unpickle( initCoordinateFile )             
                    coordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( i, j ) )
                    relaxRun = Geometry_Optimization( self.molecule, self.baseName )
                    relaxRun.Minimization()
                    Pickle( coordinateFile, self.molecule.coordinates3 ) 
            
            #...........................................................................
            with pymp.Parallel(self.nprocs) as p:
                #Pergomr the calculations for the rest of the grid
                for i in p.range ( 0, M ):
                    distance_1 = self.DINCREMENT[0] * float(i) + self.DMINIMUM[0]
                    rmodel  =  RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                    restraint_1=  RestraintMultipleDistance.WithOptions(energyModel = rmodel, distances= [ [ atom2, atom1, weight1 ],[ atom2, atom3, weight2 ] ] )
                    restraints["ReactionCoord"] = restraint_1

                    for j in range( 1,N ):
                        distance_2  = self.DINCREMENT[1] * float(j) + self.DMINIMUM[1]
                        rmodel2     = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                        restraint_2 = RestraintMultipleDistance.WithOptions(energyModel = rmodel2, distances = [ [ atom5, atom4, weight3 ],[ atom5, atom6, weight4 ] ] )
                        restraints["ReactionCoord2"] = restraint_2  
                        
                        initCoordinateFile = ""
                        if  j==0:
                            initCoordinateFile = os.path.join( self.baseName,"scan_trj" , "frame{}_{}.pkl".format( i, 0 ) )
                        elif j>0:
                            initCoordinateFile = os.path.join( self.baseName,"scan_trj" , "frame{}_{}.pkl".format( i, j-1 ) )                 
                        
                        self.molecule.coordinates3 = Unpickle( initCoordinateFile )             
                        coordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( i, j ) )
                        relaxRun = Geometry_Optimization( self.molecule, self.baseName  )
                        relaxRun.ChangeDefaultParameters( self.GeoOptPars )
                        relaxRun.Minimization()
                        Pickle( coordinateFile, self.molecule.coordinates3 )


            #Recuperate energy and reaction coordinate information   
            for i in range(M):
                for j in range(N):
                    molecule = self.molecule
                    molecule.coordinates3 = Unpickle( os.path.join( self.folder_trj , "frame{}_{}.pkl".format( i, j ) ) )
                    Energy = molecule.Energy()
                    RD1 =   molecule.coordinates3.Distance( atom1, atom2 ) - molecule.coordinates3.Distance( atom2, atom3 )   
                    RD2 =   molecule.coordinates3.Distance( atom4, atom5 ) - molecule.coordinates3.Distance( atom5, atom6 )
                    self.text += "{} {} {} {} {}\n".format( i,j,RD1,RD2,Energy-en0 )


        #-----------------------------------------------------------------------------------
        textLog = open( self.baseName+"_scan2D.log", "w" ) 
        textLog.write( self.text )
        textLog.close()
        #clean the restraints assigned to the sytem
        self.molecule.DefineRestraintModel( None )

    #--------------------------------------------------------------------------
    

#==============================================================================#
#=====================END OF CLASS FILE========================================#
#==============================================================================#