#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = pDynamo2Vismol.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

import glob, math, os, os.path

from pBabel                    import *
                                     
from pCore                     import*
                                     
from pMolecule                 import*
                              
from pMolecule.MMModel         import*
from pMolecule.NBModel         import*
   
                                     
from pMolecule.QCModel         import *
from pScientific               import *
                                     
from pScientific.Arrays        import *
                                     
from pScientific.Geometry3     import *
                                     
from pScientific.RandomNumbers import *
                                     
from pScientific.Statistics    import *
from pScientific.Symmetry      import *
                                     
from pSimulation               import *
from vModel import VismolObject

import numpy as np
#==========================================================================
def get_atom_coords_from_pdynamo_system (system, atom, frame = None):
    '''
    Function to obtain the coordinates from a given atom of a System instance.
    '''
    if frame:
        xyz = system.coordinates3[atom.index]
        ##print(atom.index, atom.label, atom.atomicNumber, atom.connections, xyz[0], xyz[1], xyz[2] )
        frame.append(float(xyz[0]))
        frame.append(float(xyz[1]))
        frame.append(float(xyz[2]))
    return frame

#==========================================================================
def get_atom_info_from_pdynamo_atom_obj (atom, sequence):
    """
    It extracts the information from the atom object, 
    belonging to pdynamo, and organizes it in the form 
    of a list that will be delivered later to build the 
    vismolObj
    
    """
    
    entityLabel = atom.parent.parent.label
    useSegmentEntityLabels = False
    if useSegmentEntityLabels:
        chainID = ""
        segID   = entityLabel[0:4]
    else:
        chainID = entityLabel[0:1]
        segID   = ""

    

    resName, resSeq, iCode = sequence.ParseLabel ( atom.parent.label, fields = 3 )
    #print(atom.index, atom.label,resName, resSeq, iCode,chainID, segID,  atom.atomicNumber, atom.connections)#, xyz[0], xyz[1], xyz[2] )
    index        = atom.index
    at_name      = atom.label
    #at_pos       = np.array([float( xyz[0]), float(xyz[1]), float(xyz[2])])
    at_pos       = np.array([float(0), float(0), float(0)])
    at_resi      = int(resSeq)
    at_resn      = resName
    at_ch        = chainID
    at_symbol    = vismolSession.vConfig.ATOM_TYPES_BY_ATOMICNUMBER[atom.atomicNumber] # at.get_symbol(at_name)
    cov_rad      = at.get_cov_rad (at_symbol)

    gridpos      = [int(at_pos[0]/gridsize), int(at_pos[1]/gridsize), int(at_pos[2]/gridsize)]
    at_occup     = 0.0
    at_bfactor   = 0.0
    at_charge    = 0.0
    
    return [index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ]

#==========================================================================
def load_pDynamo_system_from_file (filein,  gridsize = 3, vismolSession =  None, frames_only = False):
    """ Function doc """
    at  =  vismolSession.vConfig.atom_types
    
    system = ImportSystem (  filein  )
    system.BondsFromCoordinates3()
    #self.system.Summary ( )
    bonds = system.connectivity.bondIndices


    # . Get the sequence.
    sequence = getattr ( system, "sequence", None )
    if sequence is None: sequence = Sequence.FromAtoms ( system.atoms, componentLabel = "UNK.1" )
    #print (sequence)
    
    atoms           = []
    frames          = []
    
    frame = []
    for atom in system.atoms.items:
        
        frame = get_atom_coords_from_pdynamo_system (system, atom, frame)
        atoms.append(get_atom_info_from_pdynamo_atom_obj(atom, sequence))
        '''
        xyz = system.coordinates3[atom.index]
        ##print(atom.index, atom.label, atom.atomicNumber, atom.connections, xyz[0], xyz[1], xyz[2] )
        frame.append(float(xyz[0]))
        frame.append(float(xyz[1]))
        frame.append(float(xyz[2]))
        
        
        entityLabel = atom.parent.parent.label
        useSegmentEntityLabels = False
        if useSegmentEntityLabels:
            chainID = ""
            segID   = entityLabel[0:4]
        else:
            chainID = entityLabel[0:1]
            segID   = ""
    
        

        resName, resSeq, iCode = sequence.ParseLabel ( atom.parent.label, fields = 3 )
        #print(atom.index, atom.label,resName, resSeq, iCode,chainID, segID,  atom.atomicNumber, atom.connections, xyz[0], xyz[1], xyz[2] )
        index        = atom.index
        at_name      = atom.label
        at_pos       = np.array([float( xyz[0]), float(xyz[1]), float(xyz[2])])
        at_resi      = int(resSeq)
        at_resn      = resName
        at_ch        = chainID
        at_symbol    = vismolSession.vConfig.ATOM_TYPES_BY_ATOMICNUMBER[atom.atomicNumber] # at.get_symbol(at_name)
        cov_rad      = at.get_cov_rad (at_symbol)

        gridpos      = [int(at_pos[0]/gridsize), int(at_pos[1]/gridsize), int(at_pos[2]/gridsize)]
        at_occup     = 0.0
        at_bfactor   = 0.0
        at_charge    = 0.0

        
        atoms.append([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])
        '''
        
    frame = np.array(frame, dtype=np.float32)
    name = os.path.basename(filein)
    
    vismol_object  = VismolObject.VismolObject(name                           = name          ,    
                                               atoms                          = atoms         ,    
                                               vismolSession                  = vismolSession ,    
                                               bonds_pair_of_indexes          = bonds         ,    
                                               trajectory                     = [frame]       ,    
                                               auto_find_bonded_and_nonbonded = False         )

    return vismol_object


#+====================================================================================
class pDynamoSession:
    """ Class doc """
    
    def __init__ (self, vismolSession = None):
        """ Class initialiser """
        self.vismolSession  = vismolSession
        self.name           = 'pDynamo_session'
        
        self.nbModel_default         = NBModelCutOff.WithDefaults ( )
        self.fixed_color             = [0.5, 0.5, 0.5]
        self.pdynamo_distance_safety = 0.5
        
        '''self.active_id is the attribute that tells which 
        system is active for calculations in pdynamo 
        (always an integer value)'''
        self.active_id       = 0
        
        
        '''Now we can open more than one pdynamo system. 
        Each new system loaded into memory is stored in 
        the form of a dictionary, which has an int as 
        an access key.'''
        self.systems =  {
                         0:None
                         }
        
        self.systems_list = []
        self.counter      = 0
        
        
    #---------------------------------------------------------------------------------------    
    def load_a_new_pDynamo_system_from_dict (self, filesin = {}, systype = 0, name = None):
        """ Function doc """
        
        
        '''Every new system is added in the form of a 
        dict, which contains the items:'''
        psystem = {
                  'id'            : 0    ,  # access number (same as the access key in the self.systems dictionary)
                  'name'          : None ,
                  'system'        : None ,  # pdynamo system
                  
                  'vismol_object' : None ,  # Vismol object associated with the system -> is the object that will 
                                            # undergo changes when something new is requested by the interface, for example: show the QC region
                  'active'        : False, 
                  'bonds'         : None ,
                  'sequence'      : None ,
                  'qc_table'      : None ,
                  'fixed_table'   : []   ,
                   }
        
        
        if systype == 0:
            system              = ImportSystem       ( filesin['amber_prmtop'] )
            system.coordinates3 = ImportCoordinates3 ( filesin['coordinates'] )
            self.define_NBModel(_type = 1, system = system)

            
        if systype == 1:
            parameters          = CHARMMParameterFileReader.PathsToParameters (filesin['charmm_par'])
            system              = ImportSystem       ( filesin['charmm_psf'] , isXPLOR = True, parameters = parameters )
            system.coordinates3 = ImportCoordinates3 ( filesin['coordinates'] )
            self.define_NBModel(_type = 1, system = system)

        
        
        if systype == 2:
            mmModel        = MMModelOPLS.WithParameterSet ( filesin['opls_folder'] )            
            system         = ImportSystem       ( filesin['coordinates'])
            system.DefineMMModel ( mmModel )
            self.define_NBModel(_type = 1, system = system)


                
        if systype == 3:
            system = ImportSystem (filesin['coordinates'])
            system.Summary()
            print ('mmModel',system.mmModel)
            print ('qcModel',system.qcModel)
            print ('nbModel',system.nbModel)

        
        
        if name:
            pass
        else:
            name = system.label
        
        psystem['system'] =  system
        psystem['name'] =  name
        #self.name  =  name
        self.append_system_to_pdynamo_session(psystem)
        self.build_vismol_object_from_pDynamo_system (name = 'initial coordinates' )#psystem['system'].label)

    #---------------------------------------------------------------------------------------
    def append_system_to_pdynamo_session (self, psystem):
        """ Function doc """
        psystem['id']               = self.counter
        self.systems[psystem['id']] = psystem 
        
        self.systems_list.append(psystem)
        self.counter += 1
        self.active_id       = self.counter -1
    
    #---------------------------------------------------------------------------------------
    def get_bonds_from_pDynamo_system(self, safety = 0.5, id_system = False):
        self.systems[self.active_id]['system'].BondsFromCoordinates3(safety = safety)
        self.systems[self.active_id]['bonds'] = self.systems[self.active_id]['system'].connectivity.bondIndices
        return True
    
    #---------------------------------------------------------------------------------------   
    def define_NBModel (self, _type = 1 , parameters =  None, system = None):
        """ Function doc """
        
        if _type == 0:
            self.nbModel = NBModelFull.WithDefaults ( )
        
        elif _type == 1:
            self.nbModel = NBModelCutOff.WithDefaults ( )
        
        elif _type == 2:
            self.nbModel = NBModelORCA.WithDefaults ( )
        
        elif _type == 3:
            self.nbModel = NBModelDFTB.WithDefaults ( )
        
        if system:
            system.DefineNBModel ( self.nbModel )
            system.Summary ( )
        else:
            self.systems[self.active_id]['system'].DefineNBModel ( self.nbModel )
            self.systems[self.active_id]['system'].Summary ( )
        
        return True
    
    #---------------------------------------------------------------------------------------
    def get_sequence_from_pDynamo_system (self):
        """ Function doc """
        
        self.systems[self.active_id]['sequence'] = getattr ( self.systems[self.active_id]['system'], "sequence", None )
        if  self.systems[self.active_id]['sequence'] is None: 
            self.systems[self.active_id]['sequence'] = Sequence.FromAtoms ( self.systems[self.active_id]['system'].atoms, 
                                                                                              componentLabel = "UNK.1" )
        return True
    
    #---------------------------------------------------------------------------------------
    def get_atom_coords_from_pdynamo_system (self, system = None,  atom = None):

        xyz = self.systems[self.active_id]['system'].coordinates3[atom.index]
        return [float(xyz[0]),float(xyz[1]), float(xyz[2])]
    
    #---------------------------------------------------------------------------------------
    def get_atom_info_from_pdynamo_atom_obj (self, system = None, atom = None):
        """
        It extracts the information from the atom object, 
        belonging to pdynamo, and organizes it in the form 
        of a list that will be delivered later to build the 
        vismolObj
        
        """

        entityLabel = atom.parent.parent.label
        useSegmentEntityLabels = False
        if useSegmentEntityLabels:
            chainID = ""
            segID   = entityLabel[0:4]
        else:
            chainID = entityLabel[0:1]
            segID   = ""

        

        resName, resSeq, iCode = self.systems[self.active_id]['sequence'].ParseLabel ( atom.parent.label, fields = 3 )
        ##print(atom.index, atom.label,resName, resSeq, iCode,chainID, segID,  atom.atomicNumber, atom.connections)#, xyz[0], xyz[1], xyz[2] )
        
        index        = atom.index
        at_name      = atom.label
        at_resi      = int(resSeq)
        at_resn      = resName
        at_ch        = chainID
        at_symbol    = self.vismolSession.vConfig.ATOM_TYPES_BY_ATOMICNUMBER[atom.atomicNumber] # at.get_symbol(at_name)
        at_occup     = 0.0
        at_bfactor   = 0.0
        at_charge    = 0.0
        atom         = {
              'index'      : index      , 
              'name'       : at_name    , 
              'resi'       : at_resi    , 
              'resn'       : at_resn    , 
              'chain'      : at_ch      , 
              'symbol'     : at_symbol  , 
              'occupancy'  : at_occup   , 
              'bfactor'    : at_bfactor , 
              'charge'     : at_charge   
              }
        
        return atom
    
    #---------------------------------------------------------------------------------------
    def build_vismol_object_from_pDynamo_system (self                       , 
                                                 name = 'a_new_vismol_obj'  ,
                                                 system               = None,
                                                 vismol_object_active = True,
                                                 autocenter = True          ,):
        """ Function doc """
        self.get_bonds_from_pDynamo_system(safety = self.pdynamo_distance_safety)
        self.get_sequence_from_pDynamo_system()
        frames = []

        atoms = []     
        frame = []
        
        for atom in self.systems[self.active_id]['system'].atoms.items:
            xyz = self.get_atom_coords_from_pdynamo_system (atom   = atom)
            frame.append(xyz[0])
            frame.append(xyz[1])
            frame.append(xyz[2])
            
            atoms.append(self.get_atom_info_from_pdynamo_atom_obj(atom   = atom))
        

        
        frame = np.array(frame, dtype=np.float32)
        name  = os.path.basename(name)
        
        vismol_object  = VismolObject.VismolObject(name                           = name                ,    
                                                   atoms                          = atoms               ,    
                                                   vismolSession                  = self.vismolSession  ,    
                                                   bonds_pair_of_indexes          = self.systems[self.active_id]['bonds']     ,    
                                                   trajectory                     = [frame]             ,    
                                                   auto_find_bonded_and_nonbonded = False               )
        
        vismol_object.eb_sys_id = self.systems[self.active_id]['id']
        vismol_object.set_model_matrix(self.vismolSession.glwidget.vm_widget.model_mat)
        vismol_object.active = vismol_object_active
        vismol_object._get_center_of_mass(frame = 0)
        
        if self.systems[self.active_id]['system'].qcModel:
            sum_x = 0.0 
            sum_y = 0.0 
            sum_z = 0.0
            
            self.systems[self.active_id]['qc_table'] = list(self.systems[self.active_id]['system'].qcState.pureQCAtoms)
            total = 0
            
            for atom_index in self.systems[self.active_id]['qc_table']:
                atom = vismol_object.atoms[atom_index]
                
                coord = atom.coords (frame = 0)
                sum_x += coord[0]
                sum_y += coord[1]
                sum_z += coord[2]
                total+=1
                
                    
            center = np.array([sum_x / total,
                               sum_y / total, 
                               sum_z / total])
            
        else:
            center = vismol_object.mass_center

        self.systems[self.active_id]['vismol_object'] = vismol_object
        self.vismolSession.add_vismol_object_to_vismol_session (pdynamo_session  = self,
                                                                rep              = True, 
                                                                vismol_object    = vismol_object, 
                                                                autocenter       = autocenter)
        
        self.vismolSession.glwidget.vm_widget.center_on_coordinates(vismol_object, center)
        self.refresh_qc_and_fixed_representations()        
        return vismol_object
    
    #---------------------------------------------------------------------------------------
    def get_energy (self):
        """ Function doc """
        self.systems[self.active_id]['system'].Summary( )
        energy = self.systems[self.active_id]['system'].Energy( )
        return energy
    
    #---------------------------------------------------------------------------------------
    def define_free_or_fixed_atoms_from_iterable (self, fixedlist = []):
        """ Function doc """
        if fixedlist == []:
            self.systems[self.active_id]['fixed_table'] = []
            self.systems[self.active_id]['system'].freeAtoms = None
            #self.refresh_qc_and_fixed_representations()

        else:
            selection_fixed                             = Selection.FromIterable (fixedlist)
            self.systems[self.active_id]['fixed_table'] = list(selection_fixed)
            selection_free                              = selection_fixed.Complement( upperBound = len (self.systems[self.active_id]['system'].atoms ) )
        
            self.systems[self.active_id]['system'].freeAtoms = selection_free
            #self.refresh_qc_and_fixed_representations()

        self.refresh_qc_and_fixed_representations()
        return True
    
    #---------------------------------------------------------------------------------------
    def define_a_new_QCModel (self, parameters = None):
        """ Function doc """
        
        #print(parameters)
        
        electronicState = ElectronicState.WithOptions ( charge = parameters['charge'], multiplicity = parameters['multiplicity'])
        self.systems[self.active_id]['system'].electronicState = electronicState

        qcModel         = QCModelMNDO.WithOptions ( hamiltonian = parameters['method'])
        

        
        if self.systems[self.active_id]['qc_table'] :
            self.systems[self.active_id]['system'].DefineQCModel (qcModel, qcSelection = Selection.FromIterable ( self.systems[self.active_id]['qc_table']) )
            self.refresh_qc_and_fixed_representations()
            
            #print('define NBModel = ', self.nbModel)
            self.systems[self.active_id]['system'].DefineNBModel ( NBModelCutOff.WithDefaults ( ) )
        
        else:
            self.systems[self.active_id]['system'].DefineQCModel (qcModel)
            self.refresh_qc_and_fixed_representations()
    
    #---------------------------------------------------------------------------------------
    def refresh_qc_and_fixed_representations (self):
        """ Function doc >>> 
        list(molecule.qcState.boundaryAtoms) 
        list(molecule.qcState.pureQCAtoms)  self.systems[self.active_id]['system'].eh_qc_table ) )
        list(molecule.qcState.qcAtoms)
        """
        #if self.selection_fixed_table:
        #print('\n\n\nselection_fixed_table', self.systems[self.active_id]['fixed_table'])
        
        self.vismolSession.set_color_by_index(vismol_object = self.systems[self.active_id]['vismol_object'] , 
                                              indexes       = self.systems[self.active_id]['fixed_table']   , 
                                              color         = self.fixed_color)
        
        if self.systems[self.active_id]['system'].qcModel:

            self.systems[self.active_id]['qc_table'] = list(self.systems[self.active_id]['system'].qcState.pureQCAtoms)
            #boundaryAtoms     = list(self.system.qcState.boundaryAtoms)
            
            vismol_object = self.systems[self.active_id]['vismol_object']
            
            # Here we have to hide all the sticks and spheres so that there is no confusion in the representation of the QC region
            self.vismolSession.selections[self.vismolSession.current_selection].selecting_by_indexes (vismol_object = vismol_object, indexes = range(0, len(vismol_object.atoms)))
            self.vismolSession.show_or_hide_by_object (_type = 'spheres',  vobject = vismol_object, selection_table = range(0, len(vismol_object.atoms)),  show = False )
            self.vismolSession.show_or_hide_by_object (_type = 'sticks',   vobject = vismol_object, selection_table = range(0, len(vismol_object.atoms)),  show = False )
            
            self.vismolSession.selections[self.vismolSession.current_selection].unselecting_by_indexes (vismol_object = vismol_object, indexes = range(0, len(vismol_object.atoms)))
            
            self.vismolSession.selections[self.vismolSession.current_selection].selecting_by_indexes (vismol_object = vismol_object, indexes = self.systems[self.active_id]['qc_table'])
            ##print('\n\n',self.vismolSession.selections[self.vismolSession.current_selection].selected_atoms)
            
            self.vismolSession.show_or_hide_by_object (_type = 'spheres', vobject = vismol_object, selection_table = self.systems[self.active_id]['qc_table'] , show = True )
            self.vismolSession.show_or_hide_by_object (_type = 'sticks' , vobject = vismol_object, selection_table = self.systems[self.active_id]['qc_table'] , show = True )
            self.vismolSession.selections[self.vismolSession.current_selection].unselecting_by_indexes (vismol_object = vismol_object, indexes = range(0, len(vismol_object.atoms)))
            ##print('\n\n',self.vismolSession.selections[self.vismolSession.current_selection].selected_atoms)

        else:
            pass
    
    #---------------------------------------------------------------------------------------
    def import_trajectory (self, traj = None, first = 0 , last = -1, stride = 1):
        """ Function doc """
        
        traj   = '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/book/generatedFiles/cyclohexane_sdpath.ptGeo'
        frames = []
        frame  = []
        
        for atom in self.systems[self.active_id]['system'].atoms.items:
            xyz = self.get_atom_coords_from_pdynamo_system (atom   = atom)
            frame.append(xyz[0])
            frame.append(xyz[1])
            frame.append(xyz[2])
        frame = np.array(frame, dtype=np.float32)
        
        
        # . Define the trajectory.
        trajectory = ImportTrajectory ( traj, self.systems[self.active_id]['system'] )
        trajectory.ReadHeader ( )
        
        # . Loop over the frames in the trajectory.
        phi = []
        psi = []
        while trajectory.RestoreOwnerData ( ):
            frame = []
            for atom in self.systems[self.active_id]['system'].atoms.items:
                xyz = self.get_atom_coords_from_pdynamo_system (atom   = atom)
                frame.append(xyz[0])
                frame.append(xyz[1])
                frame.append(xyz[2])
            frame = np.array(frame, dtype=np.float32)
            #frames.append(frame)
            self.systems[self.active_id]['vismol_object'].frames.append(frame)

        # . Finish up.
        trajectory.ReadFooter ( )
        trajectory.Close ( )
        #return frames
    #--------------------------------------------------------------------------------------- 
    #tirar essa fũnção depois, tem os métodos de otimização no core
    def run_ConjugateGradientMinimize_SystemGeometry (self                   , 
                                                      logFrequency           , 
                                                      maximumIterations      , 
                                                      rmsGradientTolerance   , 
                                                      save_trajectory = False,
                                                      trajectory_path = None):
        """ Function doc """
        if save_trajectory:
            
            #if trajectory_path == None:
                 
            trajectory_path = '/home/fernando/Documents'
            trajectory = ExportTrajectory ('/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/book/generatedFiles/cyclohexane_sdpath.ptGeo', self.systems[self.active_id]['system'] )

            ConjugateGradientMinimize_SystemGeometry ( self.systems[self.active_id]['system']                        ,
                                                       logFrequency                       = logFrequency             ,
                                                       maximumIterations                  = maximumIterations        ,
                                                       rmsGradientTolerance               = rmsGradientTolerance     ,
                                                       trajectory                         = trajectory
                                                       )        
        
        else:
        
            ConjugateGradientMinimize_SystemGeometry ( self.systems[self.active_id]['system']                        ,
                                                       logFrequency                       = logFrequency             ,
                                                       maximumIterations                  = maximumIterations        ,
                                                       rmsGradientTolerance               = rmsGradientTolerance     )
        
        self.build_vismol_object_from_pDynamo_system (name = 'geometry optimization', autocenter = False)

        
#======================================================================================================================