#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = CoreInterface.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

import commonFunctions
import SimulationsPreset
import LogFile

#==============================================================================
class SimulationProject:
    '''
    Class to setup pDynamo simulations and methods for the system project instanced represented on vismol
    '''  
    #.-------------------------------------------------------------------------
    def __init__(self,_projectName):
        '''
        Class constructor
        '''        
        self.baseName       = projectName
        self.cSystem        = None          #current system
        self.SystemStates   = []            #list of Systems as hystory 
        self.simulation     = "Single-Point"#Name of the Simulation type
        self.simulationHys  = []            #Simulation type Hystory
        self.TheoryLevel    = "MM"          #Current state of the theory level; other types: QM and Hybrid
        self.QClist         = []
        self.BoundaryAtoms  = [] 
        self.TotalChargeMM  = 0 
        self.TotalChargeQC  = 0
        self.multiplicity   = 1 
        self.systemCoutCurr = 0 

        #Some status variable for current Sytem instance
        self.hasNBmodel     = False
        self.hasMMmodel     = False
        self.hasQCmodel     = False
        self.NBmodel        = None 
        self.QCmodel        = None
        self.MMmodel        = None
        self.ORCA           = False

        self.logfile = LogFile( _projectName + "_EasyHybrid3.log")

    #--------------------------------------------------------------------------
    def LoadSystemFromForceField(self,_topologyFile,_coordinateFile,_MDpackage):
        '''
        Class method to load the current system from topology and coordinate files from molecular dynamics packages
        '''
        if self.cSystem == None: 
            self.cSystem = ImportSystem(_topologyFile)
            self.cSystem.coordinates3 = ImportCoordinates3(_coordinateFile)
            self.systemCoutCurr += 1
        else:
            self.SystemStates.append(self.cSystem) # put the current system 
            self.cSystem = None
            self.cSystem = ImportSystem(_topologyFile)
            self.cSystem.coordinates3 = ImportCoordinates3(_coordinateFile)
            self.systemCoutCurr += 1

        #testing the MMmodel
        nbModel = NBModelCutOff.WithDefaults ( )
        self.cSystem.DefineNBModel( nbModel )

        try:
            self.cSystem.Summary()
            energy              = self.cSystem.Energy( doGradients = True )
            self.hasMMmodel     = True
            self.hasNBmodel     = True
            self.TotalChargeMM  = GetTotalCharge(self.cSystem)
        except:
            print("Problems in testing the MM model loaded from topology file!")

        self.NBmodel = self.cSystem.energyModel.nbModel
        self.MMmodel = self.cSystem.energyModel.mmModel        

    #--------------------------------------------------------------------------
    def LoadSystemFromCoordinates(self,_coordinateFile):
        '''
        Class method to Load the current system from coordinate file
        '''
        if self.cSystem == None: 
            self.cSystem = ImportSystem(_coordinateFile)
            self.systemCoutCurr += 1
        else:
            self.SystemStates.append(self.cSystem) # put the current system 
            self.cSystem = None
            self.cSystem = ImportSystem(_coordinateFile)            
            self.systemCoutCurr += 1

        nbModel = NBModelCutOff.WithDefaults ( )
        self.cSystem.DefineNBModel( nbModel )
        mmModel = MMModelOPLS.WithParameterSet ( "protein" )
        self.cSystem.DefineMMModel( mmModel )

        try:
            self.cSystem.Summary()
            energy              = self.cSystem.Energy( doGradients = True )
            self.hasMMmodel     = True
            self.hasNBmodel     = True
            self.TotalChargeMM  = GetTotalCharge(self.cSystem)
        except:
            print("Problems in testing the MM model loaded from topology file!")

        self.NBmodel = self.cSystem.energyModel.nbModel
        self.MMmodel = self.cSystem.energyModel.mmModel 
    #--------------------------------------------------------------------------
    def LoadSystemFromSavedProject(self,_pklPath):
        '''
        Class method to load the current system from a pDynamo pkl. 
        The system saved on the pkl is expected to have Non-bonded and MM models
        '''
        if self.cSystem == None: 
            self.cSystem = ImportSystem(_pklPath)
            self.systemCoutCurr += 1
        else:
            self.SystemStates.append( self.cSystem ) # put the current system 
            self.cSystem = None
            self.cSystem = ImportSystem(_pklPath)
            self.systemCoutCurr += 1

        try:
            self.cSystem.Summary()
            energy              = self.cSystem.Energy( doGradients = True )
            self.hasMMmodel     = True
            self.hasNBmodel     = True
            self.TotalChargeMM  = GetTotalCharge( self.cSystem )
        except:
            print("Problems in testing the MM model loaded from saved pkl file!")
        
        self.NBmodel = self.cSystem.energyModel.nbModel
        self.MMmodel = self.cSystem.energyModel.mmModel
        self.QCmodel = self.cSystem.energyModel.qcModel 

    #.-------------------------------------------------------------------------
    def SelectQCRegion(self,_region,_charge):
        '''
        Class method to select or change the Quantum Chemistry treated region for the system.
        If the system already has a QC energy model associated it will be kept. 

        The results of this function changes the current System instance and 
        '''

        self.TotalChargeQC = _charge        
        self.SystemStates.append( self.cSystem )
        self.cSystem.electronicState = ElectronicState ( charge = self.TotalChargeQC, multiplicity = self.multiplicity )
        self.cSystem.DefineQCModel( self.QCmodel , _region )
        self.cSystem.DefineNBModel( self.NBmodel ) # reseting the non-bonded model
        self.Summary()
        self.cSystem.Energy()


    #.-------------------------------------------------------------------------
    def SelectQCModel(self,_method,_region,_QCcharge,_QCmultiplicity):
        '''
        Class method to set a quantum chemistry Energy Model for the current system.  
        '''
        self.multiplicity   = _QCmultiplicity
        self.TotalChargeQC  = _QCcharge        
        self.SystemStates.append( self.cSystem )
        self.cSystem.electronicState = ElectronicState ( charge = self.TotalChargeQC, multiplicity = self.multiplicity )
        self.cSystem.DefineQCModel( _method , _region )
        self.cSystem.DefineNBModel( self.NBmodel ) # reseting the non-bonded model
        self.Summary()
        
        if not self.ORCA:
            self.cSystem.Energy()
        
    
    #.-------------------------------------------------------------------------
    def SetOrcaSystem(self,_model,_basis):
        '''
        Set ORCA QCMM simulation option for the current system.
        '''
        pass
    #.-------------------------------------------------------------------------
    def RunSinglePoint(self):
        '''
        Calculate the energy for the system.
        '''
        Energy = self.cSystem.Energy()
        print("Single Point Energy Calculations Done!\n")
        print(Energy)
        return(Energy)

    #.-------------------------------------------------------------------------
    def RunSimulation(self,_parameters,_simulationType):

        self.simulationHys.append(_simulationType)
        process = Simulation(self.cSystem,_simulationType)
        try:
            process.Execute(_parameters)
        except:
            Print("Simulation Type set failed to Run!")


        
   #.------------------------------------------------------------------------- 
    def Analysis(self):
        '''
        Class method to execute analysis and produce plots.
        '''
        pass

    #.-------------------------------------------------------------------------
    def RunTests(self):
        '''
        Class method to run simulation tests. 
        '''
        pass 
    
    #.-------------------------------------------------------------------------
    def SaveProject(self):
        '''
        '''
        pass


#==============================================================================

if __name__ == "__main__":
    tests = SimulationProject("projectTest")
    tests.RunTests()
