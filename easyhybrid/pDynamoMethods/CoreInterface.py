#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = CoreInterface.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

import os, glob, sys

#--------------------------------------------------------------
#Loading own libraries
from commonFunctions import *
from SimulationsPreset import Simulation
import LogFile

#--------------------------------------------------------------
#loading pDynamo Libraries
from pBabel                    import *                                     
from pCore                     import *
                                     
from pMolecule                 import *                              
from pMolecule.MMModel         import *
from pMolecule.NBModel         import *                                     
from pMolecule.QCModel         import *

from pScientific               import *                                
from pSimulation               import *

#==========================================================================

#*****************************************************************************
class SimulationProject:
    '''
    Class to setup pDynamo simulations from a remote framework, i.e. without using VisMol GUI 
    '''  
    #.-------------------------------------------------------------------------
    def __init__(self,_projectName,DEBUG=False):
        '''
        Class constructor
        DEBUG: if this paramters it set True, some extra steps in the system setting will be performe to check the things up
        '''        
        self.baseName       = _projectName
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
        self.QCRegion       = None

        self.DEBUG          = DEBUG

        #Some status variable for current Sytem instance ( For now without any great utility )
        self.hasNBmodel     = False
        self.hasMMmodel     = False
        self.hasQCmodel     = False
        self.NBmodel        = None 
        self.QCmodel        = None
        self.MMmodel        = None
        self.ORCA           = False
        self.DFTB           = False

        self.logfile = LogFile.LogFile(_projectName+"_EasyHybrid3.log")

    #===================================================================================
    def LoadSystemFromForceField(self,_topologyFile,_coordinateFile):
        '''
        Class method to load the current system from topology and coordinate files from molecular dynamics packages
        '''
        #------------------------------------------------------------
        oldSystem = Clone(self.cSystem) 
        self.logfile.separator()
        self.logfile.inputLine("None Current System Found. Creating one from topologies and coordinate files!")
        #------------------------------------------------------------
        self.cSystem = ImportSystem(_topologyFile)
        self.cSystem.coordinates3 = ImportCoordinates3(_coordinateFile)
        self.cSystem.label = self.baseName + "_#" + str(self.systemCoutCurr)
        self.systemCoutCurr += 1
        #------------------------------------------------------------
        if not oldSystem == None:
            oldSystem = copySystem(self.cSystem)
            self.SystemStates.append(oldSystem) # put the current system 
        #------------------------------------------------------------
        #testing the MMmodel
        self.nbModel = NBModelCutOff.WithDefaults ( )
        self.cSystem.DefineNBModel( self.nbModel )

        energy              = self.cSystem.Energy( doGradients = True )
        self.hasMMmodel     = True
        self.hasNBmodel     = True
        #self.TotalChargeMM  = GetTotalCharge(self.cSystem)    
        self.NBmodel = self.cSystem.nbModel
        self.MMmodel = self.cSystem.mmModel
        self.logfile.inputLine("Energy Model loaded: " + self.cSystem.energyModelLabel )
        self.logfile.inputLine("New System loaded!")


    #=====================================================================================
    def LoadSystemFromCoordinates(self,_coordinateFile):
        '''
        Class method to Load the current system from coordinate file
        #this method has import problems
        Though the use of the remote easyhtbrid without preprare topologies is unexpected
        '''
        if self.cSystem == None: 
            self.logfile.separator()
            self.logfile.inputLine("None Current System Found. Creating one from a coordinate files!")
            self.cSystem = ImportSystem(_coordinateFile)
            self.cSystem.label = self.baseName + "_#" + str(self.systemCoutCurr)
            self.systemCoutCurr += 1
        else:
            self.logfile.separator()
            self.logfile.inputLine("There is already a loaded System. Back off the current System and Creating one new from a coordinate file!")
            oldSystem = copySystem(self.cSystem)
            self.SystemStates.append(oldSystem) # put the current system 
            self.cSystem = ImportSystem(_coordinateFile) 
            self.cSystem.label = self.baseName + "_#" + str(self.systemCoutCurr)          
            self.systemCoutCurr += 1

        self.mmModel = MMModelOPLS.WithParameterSet( "protein" )
        self.cSystem.DefineMMModel( self.mmModel )
        self.nbModel = NBModelCutOff.WithDefaults ( )
        self.cSystem.DefineNBModel( self.nbModel )      

        try:
            energy              = self.cSystem.Energy( doGradients = True )
            self.hasMMmodel     = True
            self.hasNBmodel     = True
        except:
            print("Problems in testing the MM model loaded from topology file!")

        self.NBmodel = self.cSystem.nbModel
        self.MMmodel = self.cSystem.mmModel         
        self.logfile.inputLine("New System loaded!")

    #====================================================================================
    def LoadSystemFromSavedProject(self,_pklPath):
        '''
        Class method to load the current system from a pDynamo pkl. 
        The system saved on the pkl is expected to have Non-bonded and MM models
        '''
        if self.cSystem == None:
            self.logfile.separator()
            self.logfile.inputLine("None Current System Found. Creating one from a coordinate files!")
            self.cSystem = ImportSystem(_pklPath)
            self.cSystem.label = self.baseName + "_#" + str(self.systemCoutCurr)
            self.systemCoutCurr += 1
        else:
            self.logfile.separator()
            self.logfile.inputLine("There is already a loaded System. Back off the current System and Creating one new from a coordinate file!")
            oldSystem = copySystem(self.cSystem)
            self.SystemStates.append( oldSystem ) # put the current system 
            self.cSystem = ImportSystem(_pklPath)
            self.cSystem.label = self.baseName + "_#" + str(self.systemCoutCurr)
            self.systemCoutCurr += 1

        self.NBmodel = self.cSystem.nbModel
        self.MMmodel = self.cSystem.mmModel
        self.QCmodel = self.cSystem.qcModel 

        if not self.QCmodel == None:
            self.BoundaryAtoms = list(self.cSystem.qcState.boundaryAtoms)
            self.QClist        = list(self.cSystem.qcState.qcAtoms)
            self.hasQCmodel    = True

        try:
            energy              = self.cSystem.Energy( doGradients = True )
            self.hasMMmodel     = True
            self.hasNBmodel     = True
        except:
            print("Problems in testing the MM model loaded from saved pkl file!")


        self.logfile.inputLine("New System loaded!")

    #====================================================================================
    def SphericalPruning(self, _centerAtom, _radius):
        '''
        Perform a spherical pruning from a certain atom center coordinates.
        '''
        #---------------------------------------------------
        self.logfile.separator()
        self.logfile.inputLine("Starting Spherical pruning!")
        self.logfile.inputLine("BackOffing onld System!")
        #---------------------------------------------------
        oldSystem = copySystem(self.cSystem)
        self.SystemStates.append(oldSystem)
        self.systemCoutCurr += 1
        #---------------------------------------------------
        atomref = AtomSelection.FromAtomPattern( self.cSystem, _centerAtom )
        core    = AtomSelection.Within(self.cSystem,atomref,_radius)
        core2   = AtomSelection.ByComponent(self.cSystem,core)
        self.cSystem = PruneByAtom( self.cSystem , Selection(core2) )
        #---------------------------------------------------
        self.cSystem.label = self.baseName + "#{} Pruned System ".format(self.systemCoutCurr) 
        self.cSystem.DefineNBModel( self.nbModel )
        self.cSystem.Energy()

    #======================================================================================
    def SettingFixedAtoms(self,_centerAtom,_radius):
        '''
        Set the list of atoms to keep with the positions fixed through the next simulations
        '''
        #-----------------------------------------------------
        self.logfile.separator()
        self.logfile.inputLine("Setting fixed atoms!")
        self.logfile.inputLine("BackOffing onld System!")
        #-----------------------------------------------------
        oldSystem = copySystem(self.cSystem)
        self.SystemStates.append(oldSystem)
        self.systemCoutCurr += 1
        #-----------------------------------------------------
        atomref   = AtomSelection.FromAtomPattern( self.cSystem, _centerAtom )
        core      = AtomSelection.Within(self.cSystem,atomref,_radius)
        mobile    = AtomSelection.ByComponent(self.cSystem,core)        
        
        #-----------------------------------------------------
        if self.DEBUG:
            MobileSys = PruneByAtom( self.cSystem, Selection(mobile) )
            ExportSystem("MobileSystemCheck.pdb",MobileSys)
        #-----------------------------------------------------        
        self.cSystem.freeAtoms = mobile       
        self.cSystem.label = self.baseName + "#{} With Fixed Atoms ".format(self.systemCoutCurr) 
        self.cSystem.DefineNBModel( self.nbModel )
        self.cSystem.Energy()
    
    #=======================================================================================
    def ChangeQCRegion(self,_region,_QCcharge,_QCmultiplicity):
        '''
        Class method to select or change the Quantum Chemistry treated region for the system.
        If the system already has a QC energy model associated it will be kept. 

        '''
        #---------------------------------------------------------------------------------
        self.QCRegion = _region # must be already a instance o Selection class 
        #---------------------------------------------------------------------------------
        self.logfile.separator()
        self.logfile.inputLine("Changing the QC selected atoms!")
        #---------------------------------------------------------------------------------
        self.TotalChargeQC = _QCcharge  
        self.multiplicity  = _QCmultiplicity
        #---------------------------------------------------------------------------------
        #Back offing older system   
        self.logfile.inputLine("Back off the current System and Creating one!") 
        oldSystem = copySystem(self.cSystem)
        self.SystemStates.append( oldSystem )
        self.cSystem.label = self.baseName + "#{} Changed QC region".format(self.systemCoutCurr)
        self.systemCoutCurr += 1
        #--------------------------------------------------------------------------------
        self.logfile.inputLine("QC region charge: " + str(self.TotalChargeQC) )
        self.logfile.inputLine("QC region multiplicity: " + str(self.multiplicity) )        
        # Task:
        #   output the information of QC atoms in the logfile      
        self.cSystem.electronicState = ElectronicState.WithOptions ( charge = self.TotalChargeQC, multiplicity = self.multiplicity )
        self.cSystem.DefineQCModel( self.QCmodel , _region )
        self.cSystem.DefineNBModel( self.NBmodel ) # reseting the non-bonded model
        #--------------------------------------------------------------------------------
        self.BoundaryAtoms = list(self.cSystem.qcState.boundaryAtoms)
        self.QClist        = list(self.cSystem.qcState.qcAtoms)
        self.hasQCmodel    = True
        #--------------------------------------------------------------------------------
        energy = 0.0
        if not self.ORCA:
            energy = self.cSystem.Energy()
        #---------------------------------------------------------------------------------
        self.logfile.inputLine("Total Energy of the System: " + str(energy) )

    #=====================================================================================
    def SetSMOHybridModel(self,_method,_region,_QCcharge,_QCmultiplicity):
        '''
        Set a semiempirical quantum chemistry Energy Model for the current system. 
        '''
        #---------------------------------------------------------------------
        self.logfile.separator()
        self.logfile.inputLine("Defining Semiempirical method and QC atoms regions!")
        self.logfile.inputLine("\tHamiltonian {}".format(_method) )
        #_--------------------------------------------------------------------
        if not VerifyMNDOKey(_method):
            return(-1)            
        #---------------------------------------------
        # Define the QC atoms list
        atomlist = []
        for sel in _region:
            for i in range( len(sel) ):
                atomlist.append( sel[i] )
        #---------------------------------------------
        #define QC atoms selection
        converger           = DIISSCFConverger.WithOptions( energyTolerance=2.0e-4,densityTolerance=1.0e-10, maximumIterations = 500 )
        self.QCRegion       = Selection.FromIterable(atomlist)
        self.multiplicity   = _QCmultiplicity
        self.TotalChargeQC  = _QCcharge        
        #---------------------------------------------
        #Appending sytem
        self.NBmodel = self.cSystem.nbModel
        self.cSystem.nbModel = None
        oldSystem = Clone(self.cSystem)       
        self.SystemStates.append( oldSystem )        
        self.cSystem.label = self.baseName + "#{} {} Hamiltonian and QC region Set".format(self.systemCoutCurr,_method)
        self.systemCoutCurr += 1
        #Setting QC model 
        self.cSystem.electronicState = ElectronicState.WithOptions ( charge = self.TotalChargeQC, multiplicity = self.multiplicity )
        self.QCmodel = QCModelMNDO.WithOptions( hamiltonian = _method, converger=converger )
        
        #------------------------------------------------------------------------
        if self.DEBUG:
            #Export the set QC region for visual inspection        
            qcSystem = PruneByAtom(self.cSystem, self.QCRegion)
            ExportSystem(self.baseName+"_qcSystem.pdb",qcSystem)
        #------------------------------------------------------------------------
        self.cSystem.Summary()
        self.cSystem.DefineQCModel( self.QCmodel, qcSelection =self.QCRegion )
        self.cSystem.Summary()
        self.cSystem.DefineNBModel( self.NBmodel ) # reseting the non-bonded mode        
        #------------------------------------------------------------------------
        energy = self.cSystem.Energy()  
        self.logfile.inputLine("Total Energy of the System: " + str(energy) )
        #Checking status
        self.BoundaryAtoms = list(self.cSystem.qcState.boundaryAtoms)
        self.QClist        = list(self.cSystem.qcState.qcAtoms)
        self.hasQCmodel    = True   

    #=====================================================================================
    def SetOrcaSystem(self,_model,_basis,_region,_QCcharge,_QCmultiplicity):
        '''
        Set or modify the QC model to run with ORCA.
        '''
        #----------------------------------------------
        self.logfile.separator()
        self.logfile.inputLine("Defining method and QC atoms regions to run in ORCA software!")
        #seting scratch path
        _scratch = os.path.join( orcaScratchBase,self.baseName )
        if not os.path.exists(_scratch):
            os.makedirs(_scratch)
        #---------------------------------------------
        options =  "\n% output\n"
        options +=  "print [ p_mos ] 1\n"
        options +=  "print [ p_overlap ] 5\n"
        options +=  "end # output"
        #---------------------------------------------
        atomlist = []
        for sel in _region:
            for i in range( len(sel) ):
                atomlist.append( sel[i] )
        #---------------------------------------------
        #define QC atoms selection
        self.QCRegion = Selection.FromIterable(atomlist)        
        self.multiplicity   = _QCmultiplicity
        self.TotalChargeQC  = _QCcharge 
        #---------------------------------------------
        #Appending qc system
        self.cSystem.nbModel = None
        oldSystem = Clone( self.cSystem )
        self.SystemStates.append( oldSystem )
        self.cSystem.label = self.baseName + "#{} ORCA and QC region Set".format(self.systemCoutCurr)
        self.systemCoutCurr += 1
        #---------------------------------------------
        #Setting QC model
        self.cSystem.electronicState = ElectronicState.WithOptions(charge       = self.TotalChargeQC, 
                                                                   multiplicity = self.multiplicity )
        #....................................................................................................
        self.QCmodel = QCModelORCA.WithOptions( keywords        = [ _model, _basis, options ], 
                                                deleteJobFiles  = False                      ,
                                                scratch         =_scratch                    )
        #Export the set QC region for visual inspection
        #---------------------------------------------
        if self.DEBUG:
            qcSystem = PruneByAtom(self.cSystem, self.QCRegion)
            ExportSystem(self.baseName+"_qcSystem.pdb",qcSystem)        
        #---------------------------------------------
        self.NBmodel = NBModelORCA.WithDefaults()
        self.cSystem.DefineQCModel( self.QCmodel , qcSelection=self.QCRegion )
        self.cSystem.DefineNBModel(self.NBmodel)
        #---------------------------------------------
        self.BoundaryAtoms = list(self.cSystem.qcState.boundaryAtoms)
        self.QClist        = list(self.cSystem.qcState.qcAtoms)
        self.hasQCmodel    = True 


    #==========================================================================
    def SetDFTBsystem(self,_region,_QCcharge,_QCmultiplicity):
        '''
        Set or modify the QC model to run with DFTB model.
        '''
        #----------------------------------------------
        atomlist = []
        for sel in _region:
            for i in range( len(sel) ):
                atomlist.append( sel[i] )
        #---------------------------------------------
        #define QC atoms selection
        self.QCRegion = Selection.FromIterable(atomlist)
        #---------------------------------------------
        self.logfile.separator()
        self.logfile.inputLine("Defining DFTB method and QC atoms regions!")
        self.multiplicity   = _QCmultiplicity
        self.TotalChargeQC  = _QCcharge 
        #---------------------------------------------
        #Sending the system
        self.cSystem.nbModel = None
        oldSystem = Clone( self.cSystem )      
        self.SystemStates.append( oldSystem )
        self.cSystem.label = self.baseName + "#{} DFTB and QC region Set".format(self.systemCoutCurr)
        self.systemCoutCurr += 1

        self.cSystem.electronicState = ElectronicState.WithOptions ( charge = self.TotalChargeQC, multiplicity = self.multiplicity )
       
        #---------------------------------------------
        #Export the set QC region for visual inspection
        if self.DEBUG:
            qcSystem = PruneByAtom(self.cSystem, self.QCRegion)
            ExportSystem(self.baseName+"_qcSystem.pdb",qcSystem)
        #---------------------------------------------
        _scratch = os.path.join(os.getcwd(),self.baseName,"dftbjob")
        if not os.path.exists(_scratch):
            os.makedirs(_scratch)
        #--------------------------------------------------------------------
        #task adjust the parameters for customizable options
        self.QCmodel = QCModelDFTB.WithOptions ( deleteJobFiles = False   ,
                                                 randomScratch  = False   ,
                                                 scratch        = _scratch,
                                                 skfPath        = skfPath ,
                                                 useSCC         = True    )

        #-----------------------------------------------------------------------
        self.NBmodel = NBModelDFTB.WithDefaults()
        self.cSystem.DefineQCModel( self.QCmodel, qcSelection=self.QCRegion )
        self.cSystem.DefineNBModel( self.NBmodel ) # reseting the non-bonded model
        #--------------------------------------------------------------------
        energy = self.cSystem.Energy()      
        self.logfile.inputLine("Total Energy of the System: " + str(energy) )
        #--------------------------------------------------------------------
        self.BoundaryAtoms = list(self.cSystem.qcState.boundaryAtoms)
        self.QClist        = list(self.cSystem.qcState.qcAtoms)
        self.hasQCmodel    = True 

    #=========================================================================
    def RunSinglePoint(self):
        '''
        Calculate the energy for the system.
        '''
        #---------------------------------------------------------------
        self.logfile.inputLine("Single Point Energy calculation chosen!")
        #----------------------------------------------------------------
        energy = self.cSystem.Energy()
        #----------------------------------------------------------------
        print("Single Point Energy Calculations Done!\n")
        print(energy)
        #_---------------------------------------------------------------
        self.logfile.inputLine("Total Energy of the System: " + str(energy) )
        #----------------------------------------------------------------
        return(energy)

    #=========================================================================
    def RunSimulation(self,_parameters,_simulationType):
        '''
        Execute a preset simulation for the current system. 
        '''
        #----------------------------------------------------------------------        
        self.simulationHys.append(_simulationType)
        self.logfile.separator()
        self.logfile.inputLine("Setting Simulation Protocol:")
        self.logfile.inputLine( "\t{}".format(_simulationType) )
        #----------------------------------------------------------------------
        oldSystem = copySystem( self.cSystem )
        self.SystemStates.append( oldSystem )
        self.cSystem.label = self.baseName + "#{} Input for Simulation: {}".format(self.systemCoutCurr,_simulationType)
        self.systemCoutCurr += 1
        #---------------------------------------------------------------------
        bsname = os.path.join( os.getcwd(), self.baseName )
        process = Simulation(self.cSystem,_simulationType, bsname )
        process.Execute(_parameters)
              
        
    #.-------------------------------------------------------------------------
    def PrintSystems(self):
        '''
        Method to print the summary of the loaded systems 
        '''
        print("There are {} loaded systems".format( self.systemCoutCurr) )
        ctn = input("Type any key to print the Summary of the Systems, or 'N' to cancel this")
        if not ctn == "N":
            if len(self.SystemStates) > 0:
                for system in self.SystemStates:
                    system.Summary()
                    print("***************************************************")
                print("Now, printing the current system Summary:")
                self.cSystem.Summary()

            elif self.systemCoutCurr == 1:
                print("There is only the current System loaded!\n Printing its information below!")
            else: 
                print( "There are no loaded systems!")    

    #.-------------------------------------------------------------------------
    def SaveProject(self):
        '''
        The complete version of this function intends to save in pkl and another coordinate format
        the systems and trajectories worked in this simulations
        Though, in the current state only will save the current system to a pkl file
        '''
        Pickle(self.baseName+".pkl",self.cSystem)
        
    #.-------------------------------------------------------------------------
    def FinishRun(self):
        '''
        '''
        self.logfile.inputLine("Finishing simulation project using pDynamo3 methods!")
        self.logfile.close()


#==============================================================================

