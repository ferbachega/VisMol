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

    #====================================================================================    
    def Biplot(self,X,Y,type="rgrms"):
        '''
        '''
        x = np.array(X)
        y = np.array(Y)

        
        xsd   = np.std(x)
        xmean = np.mean(x)
        ysd   = np.std(y)
        ymean = np.mean(y)

        x = (x-xmean)
        y = (x-ymean)
        x = x/xsd
        y = y/ysd
        
        fig, (ax1) = plt.subplots(nrows=1)

        #----------------------------------------------------
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)
        #----------------------------------------------------
        # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
        #----------------------------------------------------
        # Setting plot type 
        pdf = ax1.scatter(x, y, c = z, s = 50, edgecolor = 'black')
        #----------------------------------------------------
        # Plot title
        #ax1.set_title('RG' + ' by ' + 'RMSD')

        # Hide right and top spines
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.yaxis.set_ticks_position('left')
        ax1.xaxis.set_ticks_position('bottom')
        
        # Set x and y limits
        xmin = x.min() 
        xmax = x.max() 
        ymin = y.min() 
        ymax = y.max()       
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)

        # Set x and y labels
        plt.ylabel("RG")
        plt.xlabel("RMSD")

        # Adding the color bar 
        colbar = plt.colorbar(pdf)
        colbar.set_label('Probability Density Function')   
        
        #printing varible stats 
        print("printing mean of RMSD: " + str( x.mean() ) )
        print("printing mean of RG: " + str( y.mean() ) )       
        plt.savefig( os.path.join( self.trajectoryNameCurr,"biplot.png") )
        plt.show()
        #.........

    #=====================================================================================
    def Analysis(self):
        '''
        Perform stuctural analysis in the production trajectory
        ''' 
        #-----------------------------------------------------------------------------       
        # . Get the atom masses and a selection for protein atoms only.
        masses  = Array.FromIterable ( [ atom.mass for atom in self.molecule.atoms ] )
        crd3    = ImportCoordinates3( os.path.join(self.trajectoryNameCurr,"frame0.pkl") )
        system  = AtomSelection.FromAtomPattern ( self.molecule, "*:*:*" )

        #------------------------------------------------------------------------------
        # . Calculate the radius of gyration.
        rg0 = crd3.RadiusOfGyration(selection = system, weights = masses)
        #------------------------------------------------------------------------------
        # . Save the starting coordinates.
        reference3 = Clone(crd3)
        #------------------------------------------------------------------------------
        # . Get the trajectory.
        trajectory = ImportTrajectory( self.trajectoryNameCurr , self.molecule )
        trajectory.ReadHeader( )
        #------------------------------------------------------------------------------
        # . Loop over the frames in the trajectory.
        rg  = []
        rms = []
        n   = []
        m   = 0 
        #-------------------------------------------------------------------------------
        while trajectory.RestoreOwnerData ( ):
            self.molecule.coordinates3.Superimpose ( reference3, selection = system, weights = masses )
            rg.append  ( self.molecule.coordinates3.RadiusOfGyration( selection = system, weights = masses ) )
            rms.append ( self.molecule.coordinates3.RootMeanSquareDeviation( reference3, selection = system, weights = masses ) )
            n.append(m)
            m+=1
        # . Set up the statistics calculations.
        rgStatistics  = Statistics(rg)
        rmsStatistics = Statistics(rms)

        #-------------------------------------------------------------------------------
        # . Save the results.        
        textLog = open( self.baseName+"_MDanalysis", "w" ) 
        #-------------------------------------------------------------------------------
        _Text = "rg0 rgMean rgSD rgMax rgMin\n"
        _Text += "{} {} {} {} {}\n".format(rg0,rgStatistics.mean,rgStatistics.standardDeviation,rgStatistics.maximum,rgStatistics.minimum )
        _Text += "rmsMean rmsSD rmsMax rmsMin\n"
        _Text += "{} {} {} {}\n".format(rmsStatistics.mean,rmsStatistics.standardDeviation,rmsStatistics.maximum,rmsStatistics.minimum )
        #-------------------------------------------------------------------------------
        _Text += "Frame RG RMS\n"
        for i in range(len(rg)):
            _Text += "{} {} {}\n".format(i,rg[i],rms[i])
        #--------------------------------------------------------------------------------
        textLog.write(_Text)
        textLog.close()
        #--------------------------------------------------------------------------------
        plt.plot(n, rg,  label = "Radius of Gyration")
        plt.xlabel("Frame")
        plt.ylabel("Radius of Gyration")
        plt.savefig( os.path.join( self.trajectoryNameCurr,"analysis_mdRG.png") )
        plt.show()

        #--------------------------------------------------------------------------------
        plt.plot(n, rms, label = "Root Mean Square")
        plt.xlabel("Frame")
        plt.ylabel("RMSD")
        plt.savefig( os.path.join( self.trajectoryNameCurr,"analysis_mdRMSD.png") )
        plt.show()        

        sns.jointplot(x=rg,y=rms,kind="kde",cmap="plasma",shade=True)
        plt.savefig( os.path.join( self.trajectoryNameCurr,"biplot.png") )

        plt.show()

        #..................
    #-------------------------------------------------------------------    
    def DistAnalysis(self,RCs):
        '''
        Perform distance analysis in the production trajectory
        '''
        #------------------------------------------------------------------------
        # . Get the trajectory.
        trajectory = ImportTrajectory( self.trajectoryNameProd , self.molecule )
        trajectory.ReadHeader()
        #------------------------------------------------------------------------
        # . Loop over the frames in the trajectory.
        rc1 = []
        rc2 = []
        energies = []
        n = []
        m=0
        #------------------------------------------------------------------------
        if len(RCs) == 2:
            while trajectory.RestoreOwnerData():
                energies.append( self.molecule.Energy(log=None) )
                if RCs[0].nAtoms == 3:
                    rc1.append( self.molecule.coordinates3.Distance(RCs[0].atoms[0], RCs[0].atoms[1]) - self.molecule.coordinates3.Distance(RCs[0].atoms[1], RCs[0].atoms[2]) )
                elif RCs[0].nAtoms == 2:
                    rc1.append( self.molecule.coordinates3.Distance(RCs[0].atoms[0], RCs[0].atoms[1]) )
                if RCs[1].nAtoms == 3:                    
                    rc2.append( self.molecule.coordinates3.Distance(RCs[1].atoms[0], RCs[1].atoms[1]) - self.molecule.coordinates3.Distance(RCs[1].atoms[1], RCs[1].atoms[2]) )
                elif RCs[1].nAtoms == 2:
                    rc2.append( self.molecule.coordinates3.Distance(RCs[1].atoms[0], RCs[1].atoms[1]) )
                n.append(m)
                m+=1
        if len(RCs) == 1:
            while trajectory.RestoreOwnerData():
                energies.append( self.molecule.Energy(log=None) )
                if RCs[0].nAtoms == 3:
                    rc1.append( self.molecule.coordinates3.Distance(RCs[0].atoms[0], RCs[0].atoms[1]) - self.molecule.coordinates3.Distance(RCs[0].atoms[1], RCs[0].atoms[2]) )
                elif RCs[0].nAtoms == 2:
                    rc1.append( self.molecule.coordinates3.Distance(RCs[0].atoms[0], RCs[0].atoms[1]) )
                n.append(m)
                m+=1

        #------------------------------------------------------------------------
        # . Save the results.        
        textLog = open( self.baseName+"_MDdistAnalysis", "w" ) 

        _Text = ""
        if len(RCs) > 1:
            _Text = "Frame RC1 RC2 Energy\n"
            for i in range(len(rc1)):
                _Text += "{} {} {} {}\n".format(i,rc1[i],rc2[i],energies[i])
        else:
            _Text = "Frame RC1 Energy\n"
            for i in range(len(rc1)):
                _Text += "{} {} {}\n".format(i,rc1[i],energies[i])

        #------------------------------------------------------------------------
        textLog.write(_Text)
        textLog.close()

        nenergies = []
        for i in range(len(energies) ):
            nenergies.append(energies[i] - energies[0])

        #---------------------------------------------

        plt.plot(n, nenergies)
        plt.savefig(self.trajectoryNameCurr+"_MDenergy.png")
        plt.show()
        
        plt.plot(n, rc1)
        if len(RCs) ==2:
            plt.plot(n, rc2, label = "Distance #2")
        plt.savefig(self.trajectoryNameCurr+"_MDdistAnalysis.png")
        plt.show()

        if len(RCs) == 2:
            sns.jointplot(x=rc1,y=rc2,kind="kde",cmap="plasma",shade=True)
            plt.savefig( os.path.join( self.trajectoryNameCurr,"distanceBiplot.png") )
            plt.show()
        #................................................        

#===============================================================================#
#. End of class MD =============================================================#
#===============================================================================#
       
    
    