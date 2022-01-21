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


#---------------------------------------

#==============================================================================
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
        if not os.path.exists(_baseFolder):
            os.makedirs(_baseFolder)
        self.trajectoryNameEqui     = os.path.join(_baseFolder,"equilibration.ptGeo")
        self.trajectoryNameProd     = os.path.join(_baseFolder,"production.ptGeo")
        self.trajectoryNameSoft     = os.path.join(_baseFolder,"soft_constr.ptRes")
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
        self.logname                = ""
        self.softConstraint         = False # boolean flags signlizing whether the system has soft constraints or not              
        #Parameters to be used internally for information
        self.state                  = "non-setted"  # string sinalizando o status da instancia dessa classe. 
        self.equilibrated           = False         # Signlizes the whether the system passed through a succesfull equilibrtion simulation 
        self.isQuantum              = False         # Boolean flag signilizing if some extent of the system is treated with a quatum energy model
        self.equi_time              = self.equiNsteps * self.timeStep            # default time for equilibration in ps
        self.prod_time              = self.prodNsteps * self.timeStep # production time where the data must to be collected
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
    
    #.---------------------------------------    
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
        if 'save_frequency'             in _parameters: 
            self.samplingFactor         = _parameters['saveFreq']
        if 'log_frequency'              in _parameters:
            self.logFrequency           = _parameters['log_frequency']
        if 'temperature_scale_option'   in _parameters:
            self.temperatureScaleOption = _parameters['temperature_scale_option']

    #.---------------------------------------    
    def HeatingSystem(self,_nsteps):
        '''
        Run a Velocity Verlet molecular dynamics simulation to gradually 
        make the system reach certain temperature. 
        '''
        
        VelocityVerletDynamics_SystemGeometry(  self.molecule                                        ,
                                                logFrequency              = self.logFreq             ,
                                                normalDeviateGenerator    = self.RNG                 ,
                                                steps                     = _nsteps                  ,
                                                timeStep                  = self.timeStep            ,
                                                temperatureScaleFrequency = self.temperatureScaleFreq,
                                                temperatureScaleOption    = self.temperatureScaleOption,
                                                temperatureStart          = self.startTemperature    ,
                                                temperatureStop           = self.temperature         )
    
        self.state = "Heated"
    #.---------------------------------------
    def RunEquilibration(self,_equiSteps):
        '''
        Run a molecular dynamics simulation for equilibration of the system
        '''

        self.nsteps             = _equiSteps
        self.trajectoryNameCurr = self.trajectoryNameEqui

        if self.algorithm == "Verlet":
            self.runVerlet()
        elif self.algorithm == "LeapFrog":
            self.runLeapFrog()
        elif self.algorithm == "Langevin":
            self.runLangevin()

            
    #.---------------------------------------    
    def RunProduction(self,_prodSteps,):
        '''
        Run a molecular dynamics simulation for data collection.
        '''
        self.nsteps             = _prodSteps
        self.trajectoryNameCurr = self.trajectoryNameProd         
        
        if self.algorithm == "Verlet":
            self.runVerlet()
        elif self.algorithm == "LeapFrog":
            self.runLeapFrog()
        elif self.algorithm == "Langevin":
            self.runLangevin()

        if self.outputDCD:
            Duplicate(self.trajectoryNameCurr,self.trajectoryNameCurr+".dcd",self.molecule)
  
   
    #.---------------------------------------
    def RunProductionRestricted(self,_equiSteps,_prodSteps,_samplingFactor):
        '''
        Run a simulation with the system having soft constrains defined.
        Designed for Umbrella Sampling and Steered Molecular Dynamics routines.
        '''
        self.nsteps             = _prodSteps
        self.samplingFactor     = _samplingFactor
        self.softConstraint     = True
        
        self.RunEquilibration(_equiSteps)
        self.RunProduction(_prodSteps)

        if self.outputDCD:
            Duplicate(self.trajectoryNameCurr,self.trajectoryNameCurr+".dcd",self.molecule)
    
    #.---------------------------------------
    def CheckSimulation(self,_type):
        '''
        Class method to check whether the molecular dynamics ended well. 
        In the case of the equilibration runs, check whether the simulations reached satisfactorly conditions. 

        Tasks: 
            1. Develop algorithm that decides if the simulation go well
            2. attach plotting functionalities for visualization 
            3. write the method
        '''
        pass   

    #.---------------------------------------
    def runVerlet(self):
        '''
        Execute velocity verlet molecular dynamics from pDynamo methods. 
        '''
        trajectory      = ExportTrajectory( self.trajectoryNameCurr, self.molecule )         
        trajectory_list = list

        if self.softConstraint:
            trajSoft = ExportTrajectory(self.trajectoryNameSoft, self.molecule)
            trajectory_list = [ (trajectory, self.samplingFactor ), (trajSoft, 1) ]
        else:
            trajectory_list = [ (trajectory, self.samplingFactor ) ]
        
        VelocityVerletDynamics_SystemGeometry(self.molecule                             ,
                                logFrequency                = self.logFreq              ,
                                normalDeviateGenerator      = self.RNG                  ,
                                steps                       = self.nsteps               ,
                                timeStep                    = self.timeStep             ,
                                temperatureScaleFrequency   = self.temperatureScaleFreq ,
                                temperatureScaleOption      = "constant"                ,
                                trajectories                = trajectory_list           ,
                                temperatureStart            = self.temperature          )

    #.---------------------------------------
    def runLeapFrog(self):
        '''
        Execute Leap Frog molecular dynamics from pDynamo methods.
        '''

        trajectory  = ExportTrajectory(self.trajectoryNameCurr, self.molecule)       
        trajectory_list = list

        if self.softConstraint:
            trajSoft = ExportTrajectory(self.trajectoryNameSoft, self.molecule)
            trajectory_list = [ ( trajectory, self.samplingFactor ), ( trajSoft, 1 ) ]
        else:
            trajectory_list = [ ( trajectory, self.samplingFactor ) ]

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

    
    #.---------------------------------------    
    def runLangevin(self):
        '''
        Execute Langevin molecular dynamics from pDynamo methods.
        '''
        trajectory  = ExportTrajectory(self.trajectoryNameCurr, self.molecule)
        trajectory_list = list

        if self.softConstraint:
            trajSoft = ExportTrajectory(self.trajectoryNameSoft, self.molecule)
            trajectory_list = [ ( trajectory, self.samplingFactor ), ( trajSoft, 1 ) ]
        else:
            trajectory_list = [ ( trajectory, self.samplingFactor ) ]

        LangevinDynamics_SystemGeometry ( self.molecule                         ,
                                          collisionFrequency     = self.collFreq,
                                          logFrequency           = self.logFreq ,
                                          normalDeviateGenerator = self.RNG     ,
                                          steps                  = self.nsteps  ,
                                          temperature            =   300.0      ,
                                          timeStep               =   0.001      ,
                                          trajectories = trajectory_list        )

    #_--------------------------------------------------------
    def Analysis(self):
        '''
        '''        
        # . Get the atom masses and a selection for protein atoms only.
        masses  = Array.FromIterable ( [ atom.mass for atom in self.molecule.atoms ] )
        crd3    = ImportCoordinates3( os.path.join(self.trajectoryNameProd,"frame1.pkl") )
        system  = AtomSelection.FromAtomPattern ( self.molecule, "A:*:*" )

        # . Calculate the radius of gyration.
        rg0 = crd3.RadiusOfGyration ( selection = system, weights = masses )

        # . Save the starting coordinates.
        reference3 = Clone ( crd3 )

        # . Get the trajectory.
        trajectory = ImportTrajectory (  self.trajectoryNameProd , self.molecule )
        trajectory.ReadHeader ( )

        # . Loop over the frames in the trajectory.
        rg  = []
        rms = []
        n   = []
        m   = 0 
        while trajectory.RestoreOwnerData ( ):
            self.molecule.coordinates3.Superimpose ( reference3, selection = system, weights = masses )
            rg.append  ( self.molecule.coordinates3.RadiusOfGyration        (             selection = system, weights = masses ) )
            rms.append ( self.molecule.coordinates3.RootMeanSquareDeviation ( reference3, selection = system, weights = masses ) )
            n.append(m)
            m+=1
        # . Set up the statistics calculations.
        rgStatistics  = Statistics ( rg  )
        rmsStatistics = Statistics ( rms )

        # . Save the results.        
        textLog = open( self.baseName+"_MDanalysis", "w" ) 

        _Text = "rg0 rgMean rgSD rgMax rgMin\n"
        _Text += "{} {} {} {} {}".format(rg0,rgStatistics.mean,rgStatistics.standardDeviation,rgStatistics.maximum,rgStatistics.minimum )
        _Text += "rms0 rmsMean rmsSD rmsMax rmsMin\n"
        _Text += "{} {} {} {} {}".format(rms0,rmsStatistics.mean,rmsStatistics.standardDeviation,rmsStatistics.maximum,rmsStatistics.minimum )

        _Text += "RG RMS\n"
        for i in range(len(rg)):
            _Text += "{} {}\n".format(rg,rms)

        textLog.write()
        textLog.close()

        plt.plot(n, rg,  label = "Radius of Gyration")
        plt.plot(n, rms, label = "Root Mean Square")
        plt.show()
        plt.savefig(self.baseName+"_mdAnalysis.png")

    
    #-------------------------------------------------------------------    
    def DistAnalysis(self,rcs,mdis):
        '''
        '''
        atom01 = rcs[0][0]
        atom02 = rcs[0][1]
        atom03 = 0 
        
        if mdis:
            atom03 = rcs[0][2]

        atom11 = 0
        atom12 = 0
        atom13 = 0 

        if mdis:
            atom03 = rcs[0][2]       

        if len(rcs) > 1:
            atom11 = rcs[1][0]
            atom12 = rcs[1][1]

            if mdis:
                atom13 = rcs[1][2]

        # . Get the trajectory.
        trajectory = ImportTrajectory (  self.trajectoryNameProd , self.molecule )
        trajectory.ReadHeader ( )

        # . Loop over the frames in the trajectory.
        rc1 = []
        rc2 = []
        energies = []
        n = []

        if mdis and len(rcs) > 1:
            while trajectory.RestoreOwnerData ( ):
                self.energies.append( self.molecule.Energy() )
                rc1.append( self.molecule.coordinates3.Distance ( atom01, atom02 ) - self.molecule.coordinates3.Distance ( atom02, atom03 ) )
                rc2.append( self.molecule.coordinates3.Distance ( atom11, atom12 ) - self.molecule.coordinates3.Distance ( atom12, atom13 ) )
                n.append(m)
                m+=1
        # . Set up the statistics calculations.

        # . Save the results.        
        textLog = open( self.baseName+"_MDdistAnalysis", "w" ) 

        _Text = "rg0 rgMean rgSD rgMax rgMin\n"
        _Text += "{} {} {} {} {}".format(rg0,rgStatistics.mean,rgStatistics.standardDeviation,rgStatistics.maximum,rgStatistics.minimum )
        _Text += "rms0 rmsMean rmsSD rmsMax rmsMin\n"
        _Text += "{} {} {} {} {}".format(rms0,rmsStatistics.mean,rmsStatistics.standardDeviation,rmsStatistics.maximum,rmsStatistics.minimum )

        _Text += "RG RMS\n"
        for i in range(len(rg)):
            _Text += "{} {}\n".format(rg,rms)

        textLog.write()
        textLog.close()

        plt.plot(n, rg,  label = "Radius of Gyration")
        plt.plot(n, rms, label = "Root Mean Square")
        plt.show()
        plt.savefig(self.baseName+"_mdAnalysis.png")

        pass

#==============================================
#. End of class MD
#==============================================
       
    
    