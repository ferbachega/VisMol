#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#FILE = RelaxedScan.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#==============================================================================

import pymp
import pymp.shared
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
        self.text = ""
            
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
        if "forceConst" in _parameters:
            self.forceC = _parameters["forceConst"]

        #-----------------------------------------------------------
        self.GeoOptPars =   { 
                            "method": self.optmizer           ,\
                            "logFrequency": self.logFreq      ,\
                            "saveTraj": "false"               ,\
                            "maxIterations":self.maxIt        ,\
                            "rmsGradient": self.rmsGT
                            }
    #=============================================================================================
    def SetReactionCoord(self,_atoms,_dincre,_massConstraint):
        '''
        Set reaction coordinate, determining initial parameters from the atoms information
        '''
        #------------------------------------------------------------
        ndim = self.nDim # temp var to hold the index of the curren dim
        self.nDim += 1
        self.atoms.append(_atoms)
        self.DINCREMENT[ndim]   = _dincre
        self.sigma_a1_a3[ndim]  =  1.0
        self.sigma_a3_a1[ndim]  = -1.0
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
        distance = 0.0
        en0 = self.molecule.Energy()

        #----------------------------------------------------------------------------------------
        if self.multipleDistance[0]:
            for i in range(_nsteps):     
                distance = self.DMINIMUM[0] + self.DINCREMENT[0] * float( i )
                rmodel = RestraintEnergyModel.Harmonic( distance, self.forceC )
                restraint = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom2, atom1, weight1 ], [ atom2, atom3, weight2 ] ] )
                restraints['ReactionCoord'] =  restraint            

                relaxRun = GeometrySearcher(self.molecule, self.baseName  )
                relaxRun.ChangeDefaultParameters(self.GeoOptPars)
                relaxRun.Minimization(self.optmizer)
                
                self.energies.append( self.molecule.Energy() - en0 )
                self.reactionCoordinate1.append( self.molecule.coordinates3.Distance( atom1 , atom2  ) - self.molecule.coordinates3.Distance( atom2, atom3  ) ) 
                Pickle( os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}.pkl".format(i) ), self.molecule.coordinates3 )
                self.text += "{} {} {} \n".format( i, self.reactionCoordinate1[i],self.energies[i]) 

        #..........................................................................................
        else:
             for i in range(_nsteps):       
                distance = self.DMINIMUM[0] + self.DINCREMENT[0] * float( i )
                rmodel = RestraintEnergyModel.Harmonic( distance, self.forceC )
                restraint = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom1, point2= atom2 )
                restraints['ReactionCoord'] =  restraint            
                #--------------------------------------------------------------------
                relaxRun = GeometrySearcher(self.molecule,"test")
                relaxRun.ChangeDefaultParameters(self.GeoOptPars)
                relaxRun.Minimization(self.optmizer)
                #--------------------------------------------------------------------
                self.energies.append( self.molecule.Energy() -en0 )
                self.reactionCoordinate1.append( self.molecule.coordinates3.Distance( atom1 , atom2  ) )  
                Pickle( os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}.pkl".format(i) ), self.molecule.coordinates3 ) 
                self.text += "{} {} {} \n".format( i, self.reactionCoordinate1[i],self.energies[i]) 
        
        
    
    #===================================================================================================
    def Run2DScan(self,_nsteps_x,_nsteps_y):
        '''
        '''
        #------------------------------------------------------
        self.text += "x y RC1 RC2 Energy\n" 
        #------------------------------------------------------
        self.restraints = RestraintModel( )
        self.molecule.DefineRestraintModel( self.restraints )
        #------------------------------------------------------
        distance_1 = 0.0 
        distance_2 = 0.0
        #------------------------------------------------------
        X = _nsteps_x
        Y = _nsteps_y
        self.nsteps[0] = X       
        self.nsteps[1] = Y  
        self.energiesMatrix = np.zeros( (X,Y),dtype=float)    
        #-----------------------------------------------------------------------------------------------
        #Define the origin point of the relaxed surface scan, AKA the 0,0 point
        coordinateFile = os.path.join( self.baseName ,"ScanTraj.ptGeo","frame{}_{}.pkl".format( 0, 0 ) )
        relaxRun = GeometrySearcher( self.molecule, self.baseName )
        relaxRun.ChangeDefaultParameters( self.GeoOptPars )
        relaxRun.Minimization(self.optmizer)
        #-----------------------------------------------------------------------------------------------
        self.en0 = self.molecule.Energy(log=None)
        Pickle( coordinateFile, self.molecule.coordinates3 )

        if self.multipleDistance[0] and self.multipleDistance[1]:
            self.Run2DScanMultipleDistance(X,Y)            
        elif self.multipleDistance[0] and self.multipleDistance[1] == False:            
            self.Run2DMixedDistance(X,Y)
        else:
            self.Run2DSimpleDistance(X,Y)

    
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
                coordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( i, 0 ) )
                relaxRun = GeometrySearcher( self.molecule, self.baseName )
                relaxRun.ChangeDefaultParameters(self.GeoOptPars)
                relaxRun.Minimization(self.optmizer)            
                #-----------------------------------------------------------------------------------
                Pickle( coordinateFile, self.molecule.coordinates3 )
                #....................................................
            #-------------------------------------------------------------------------------------------
        with pymp.Parallel(self.nprocs) as p:
            #Pergomr the calculations for the rest of the grid
            for i in p.range ( 0, X ):

                distance_1 = self.DMINIMUM[0] + ( self.DINCREMENT[0] * float(i) )
                rmodel     = RestraintEnergyModel.Harmonic( distance_1, self.forceC )
                restraint  = RestraintDistance.WithOptions(energyModel =rmodel, point1=atom1, point2=atom2  )
                restraints["RC1"] = restraint
                
                for j in range( 1, Y ):

                    distance_2  = self.DMINIMUM[1] + ( self.DINCREMENT[1] * float(j) )
                    rmodel      = RestraintEnergyModel.Harmonic( distance_2, self.forceC )
                    restraint   = RestraintDistance.WithOptions(energyModel = rmodel, point1=atom3, point2=atom4  )
                    restraints["RC2"] = restraint                    
                    
                    initCoordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo" , "frame{}_{}.pkl".format( i, j-1 ) )     
                    self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile, log=None )             
                    coordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( i, j ) )
                    relaxRun = GeometrySearcher( self.molecule, self.baseName  )
                    relaxRun.ChangeDefaultParameters( self.GeoOptPars )
                    try:
                       relaxRun.Minimization(self.optmizer)
                    except:
                        print("Possible Convergence Problems!")
                    #-----------------------------------------------------------------------------------
                    Pickle( coordinateFile, self.molecule.coordinates3 )
                    #-----------------------------------------------------------------------------------
        self.molecule.DefineRestraintModel(None)
        #-------------------------------------------------------------------------------------------------------
        for i in range(X):
            for j in range(Y):
                coordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( i, j ) )
                self.molecule.coordinates3 = ImportCoordinates3(coordinateFile,log=None)
                self.energies.append( self.molecule.Energy(log=None)-self.en0  )
                self.reactionCoordinate1.append( self.molecule.coordinates3.Distance( atom1, atom2 ) )
                self.reactionCoordinate2.append( self.molecule.coordinates3.Distance( atom3, atom4 ) )
                self.text += "{} {} {} {} {}\n".format( i,j,self.reactionCoordinate1[-1],self.reactionCoordinate2[-1],self.energies[-1])
                self.energiesMatrix[j][i] = self.energies[-1]


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

        with pymp.Parallel(self.nprocs) as p:
            for i in p.range ( 1, X ):  
            #.----              
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
                coordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( i, 0 ) )
                relaxRun = GeometrySearcher( self.molecule, self.baseName )
                relaxRun.ChangeDefaultParameters( self.GeoOptPars )
                relaxRun.Minimization(self.optmizer)                   
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
                    initCoordinateFile = ""
                    if  i==0:
                        initCoordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo" , "frame{}_{}.pkl".format( 0, j-1 ) )
                    elif i>0:
                        initCoordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo" , "frame{}_{}.pkl".format( i, j-1 ) )                 
                    #-----------------------------------------------------------------------------------
                    self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile, log=None )             
                    coordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( i, j ) )
                    relaxRun = GeometrySearcher( self.molecule, self.baseName  )
                    relaxRun.ChangeDefaultParameters( self.GeoOptPars )
                    relaxRun.Minimization(self.optmizer)
                    #-----------------------------------------------------------------------------------
                    Pickle( coordinateFile, self.molecule.coordinates3 )
                    #...................................................
        self.molecule.DefineRestraintModel(None)
        #---------------------------------------------------------------------------------------------------------
        for i in range(X):
            for j in range(Y):
                coordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( i, j ) )
                self.molecule.coordinates3 = ImportCoordinates3(coordinateFile,log=None)
                self.energies.append( self.molecule.Energy(log=None )-self.en0 )
                self.reactionCoordinate1.append( self.molecule.coordinates3.Distance( atom1, atom2 ) - self.molecule.coordinates3.Distance( atom2, atom3 ) )
                self.reactionCoordinate2.append( self.molecule.coordinates3.Distance( atom4, atom5 )  )
                self.text += "{} {} {} {} {}\n".format( i,j,self.reactionCoordinate1[-1],self.reactionCoordinate2[-1],self.energies[-1])
                self.energiesMatrix[j][i] = self.energies[-1]
    #=======================================================
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
                relaxRun = GeometrySearcher( self.molecule, self.baseName )
                relaxRun.ChangeDefaultParameters(self.GeoOptPars)
                relaxRun.Minimization(self.optmizer)
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
                    #---- ----------------------------------------------------------------------------
                    self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile, log = None )             
                    relaxRun = GeometrySearcher( self.molecule, self.baseName  )
                    relaxRun.ChangeDefaultParameters(self.GeoOptPars)
                    try:
                        relaxRun.Minimization( self.optmizer )
                    except:
                        print("Possible Convergence Problems")
                    #-----------------------------------------------------------------------------------
                    coordinateFile = os.path.join( self.baseName, "ScanTraj.ptGeo", "frame"+str(i)+"_"+str(j)+".pkl" )
                    Pickle( coordinateFile, self.molecule.coordinates3 )
                    
        #--------------------------------------                
        self.molecule.DefineRestraintModel(None)
        #-----------------------------------------------------------------------------------
        for i in range(X):
            for j in range(Y):
                coordinateFile = os.path.join( self.baseName,"ScanTraj.ptGeo", "frame{}_{}.pkl".format( i, j ) )
                self.molecule.coordinates3 = ImportCoordinates3(coordinateFile,log=None)
                self.energies.append( self.molecule.Energy(log=None)- self.en0  )
                self.reactionCoordinate1.append( self.molecule.coordinates3.Distance( atom1, atom2 ) - self.molecule.coordinates3.Distance( atom2, atom3 ) )
                self.reactionCoordinate2.append( self.molecule.coordinates3.Distance( atom4, atom5 ) - self.molecule.coordinates3.Distance( atom5, atom6 ) )
                self.text += "{} {} {} {} {}\n".format( i,j,self.reactionCoordinate1[-1],self.reactionCoordinate2[-1],self.energies[-1])
                self.energiesMatrix[j][i] = self.energies[-1]

    #======================================================================================
    def Plot1D(self):
        '''
        Make one-dimensional energy plot
        '''        
        plt.plot(self.reactionCoordinate1,self.energies)
        plt.xlabel('Reaction Coordinate (Ang)')
        plt.ylabel('Potential Energy')
        plt.savefig(self.baseName+"_1Denergy.png")
        plt.show()

    #======================================================================================
    def PlotContour(self):
        '''
        Make 2D energy plot
        '''
        x = range(self.nsteps[0])
        y = range(self.nsteps[1])        
        X = np.linspace(self.reactionCoordinate1[0],self.reactionCoordinate1[-1],self.nsteps[0])
        Y = np.linspace(self.reactionCoordinate2[0],self.reactionCoordinate2[-1],self.nsteps[1])
        
        z = self.energiesMatrix

        fig, (ax0) = plt.subplots(nrows=1)
        #setting plot parameters
        vmin=z.min()
        vmax=z.max()
        levels = MaxNLocator(nbins=100).tick_values( z.min(), z.max() )
        cmap = plt.get_cmap("jet")

        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        norm= colors.PowerNorm(gamma=1./2.)
        norm= colors.Normalize(vmin=vmin, vmax=vmax)

        im = ax0.pcolormesh(X,Y,z, cmap=cmap, norm=norm, shading = "gouraud")
        am = ax0.contour(X,Y,z,10, colors='k')
        
        ax0.clabel(am, inline=1, fontsize=8, fmt='%1.1f',colors="k")
        
        cbar = fig.colorbar(im, ax=ax0)
        cbar.ax.tick_params()

        FontSize = 14
        # Set the tick labels font
        axis_font = {'fontname':'Michroma', 'size':14}
        for tick in (ax0.xaxis.get_major_ticks()):
            tick.label.set_fontname('Arial')
            tick.label.set_fontsize(FontSize)

        for tick in (ax0.yaxis.get_major_ticks()):
            tick.label.set_fontname('Dejavu')
            tick.label.set_fontsize(FontSize) 
        

        ax0.set_xlabel("Reaction Coordinate #1", **axis_font)
        ax0.set_ylabel("Reaction Coordinate #2", **axis_font)

        fig.tight_layout()
        plt.savefig(self.baseName+"_PES.png",dpi=1000)
        #plt.show()

    #=======================================================================================
    def Finalize(self):
        '''
        Writing logs, making plots and saving trajectories
        '''
        #-----------------------------------------------------------------
        self.molecule.DefineRestraintModel( None )
        if self.nDim == 1:
            #..................................................
            trajNameDCD = self.baseName + ".dcd"
            trajName = os.path.join( self.baseName, "ScanTraj.ptGeo" )
            Duplicate( trajName, trajNameDCD, self.molecule )
            #..................................................
            textLog = open( self.baseName+"_scan1D.log", "w" ) 
            textLog.write(self.text)
            textLog.close() 
            #..................................................
            self.Plot1D() 
        #-----------------------------------------------------------------   
        elif self.nDim == 2:
            textLog = open(self.baseName+"_sca2D.log", "w") 
            textLog.write(self.text)
            textLog.close()
            #...................................................
            self.PlotContour()

#==============================================================================#
#=====================END OF CLASS FILE========================================#
#==============================================================================#
