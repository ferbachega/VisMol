#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Lembrar de colocar uma header nesse arquivo

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

import glob, math, os, os.path, sys
import numpy as np
VISMOL_HOME = os.environ.get('VISMOL_HOME')
#path fo the core python files on your machine
sys.path.append(os.path.join(VISMOL_HOME,"easyhybrid/pDynamoMethods") )
sys.path.append(os.path.join(VISMOL_HOME,"easyhybrid/gui"))

#---------------------------------------
from pBabel                    import*                                     
from pCore                     import*  
#---------------------------------------
from pMolecule                 import*                              
from pMolecule.MMModel         import*
from pMolecule.NBModel         import*                                     
from pMolecule.QCModel         import*
#---------------------------------------
from pScientific               import*                                     
from pScientific.Arrays        import*                                     
from pScientific.Geometry3     import*                                     
from pScientific.RandomNumbers import*                                     
from pScientific.Statistics    import*
from pScientific.Symmetry      import*
#---------------------------------------                              
from pSimulation               import*
#---------------------------------------
#import our core lib
from SimulationsPreset import Simulation 
#---------------------------------------
from vModel import VismolObject
from vModel.MolecularProperties import ATOM_TYPES_BY_ATOMICNUMBER
from vModel.MolecularProperties import COLOR_PALETTE

from easyhybrid.gui import *
from easyhybrid.gui.PES_analisys_window  import  PotentialEnergyAnalysisWindow 
from easyhybrid.gui.PES_analisys_window  import  parse_2D_scan_logfile 


HOME = os.environ.get('HOME')

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
    at_symbol    = ATOM_TYPES_BY_ATOMICNUMBER[atom.atomicNumber] # at.get_symbol(at_name)
    cov_rad      = at.get_cov_rad (at_symbol)

    gridpos      = [int(at_pos[0]/gridsize), int(at_pos[1]/gridsize), int(at_pos[2]/gridsize)]
    at_occup     = 0.0
    at_bfactor   = 0.0
    at_charge    = 0.0
    
    return [index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ]

#==========================================================================
def load_pDynamo_system_from_file (filein,  gridsize = 3, vm_session =  None, frames_only = False):
    """ Function doc """
    at  =  vm_session.vConfig.atom_types
    
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
        at_symbol    = vm_session.vConfig.ATOM_TYPES_BY_ATOMICNUMBER[atom.atomicNumber] # at.get_symbol(at_name)
        cov_rad      = at.get_cov_rad (at_symbol)

        gridpos      = [int(at_pos[0]/gridsize), int(at_pos[1]/gridsize), int(at_pos[2]/gridsize)]
        at_occup     = 0.0
        at_bfactor   = 0.0
        at_charge    = 0.0

        
        atoms.append([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])
        '''
        
    frame = np.array(frame, dtype=np.float32)
    name  = os.path.basename(filein)
    
    vismol_object  = VismolObject.VismolObject(name                           = name          ,    
                                               atoms                          = atoms         ,    
                                               vm_session                     = vm_session ,    
                                               bonds_pair_of_indexes          = bonds         ,    
                                               trajectory                     = [frame]       ,    
                                               auto_find_bonded_and_nonbonded = False         )

    return vismol_object


#+====================================================================================
class pDynamoSession:
    """ Class doc """
    
    def __init__ (self, vm_session = None):
        """ Class initialiser """
        self.vm_session  = vm_session
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
        
        #self.systems_list = []
        self.counter      = 0
        self.color_palette_counter = 0
        #self.sel_name_counter = 0
    
    
    def restart_pdynamo2vismol_session (self):
        """ Function doc """
        self.name           = 'pDynamo_session'
        #self.nbModel_default         = NBModelCutOff.WithDefaults ( )
        #self.fixed_color             = [0.5, 0.5, 0.5]
        #self.pdynamo_distance_safety = 0.5
        
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
        
        #self.systems_list = []
        self.counter      = 0
        self.color_palette_counter = 0
        #self.sel_name_counter = 0
    
    
    def export_system (self, sys_id, filename, folder, _format):
        """ Function doc """
        system = self.systems[sys_id]['system']
        ExportSystem ( os.path.join ( folder, filename+'.'+_format), system )

    def generate_pSystem_dictionary (self, system, working_folder = None ):
        """ Function doc """
        
        if working_folder:
            pass
        else:
            working_folder = HOME
        
        psystem = {
                  'id'            : 0              ,  # access number (same as the access key in the self.systems dictionary)
                  'name'          : None           ,
                  'system'        : system         ,  # pdynamo system
                                                   
                  'vismol_object' : None           ,  # Vismol object associated with the system -> is the object that will 
                                                      # undergo changes when something new is requested by the interface, for example: show the QC region
                  'active'        : False          , 
                  'bonds'         : None           ,
                  'sequence'      : None           ,
                  'qc_table'      : None           ,
                  'color_palette' : None           , # will be replaced by a dict
                  'fixed_table'   : []             ,
                  'selections'    : {}             ,
                  'logfile_data'  : {}             , # <--- vobject_id : [['logfilename', data], ['logfilename', data]]
                  'working_folder': working_folder , 
                   }
        
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
                  'color_palette' : None , # will be replaced by a dict
                  'fixed_table'   : []   ,
                  'selections'    : {}   ,
                  'vismol_objects': {}   ,
                  'logfile_data'  : {}   ,
                  'working_folder': HOME , 
                   }
        
        system = None 
        if systype == 0:
            system              = ImportSystem       ( filesin['amber_prmtop'] )
            system.coordinates3 = ImportCoordinates3 ( filesin['coordinates'] )
            self.define_NBModel(_type = 1, system = system)                      
        elif systype == 1:
            parameters          = CHARMMParameterFileReader.PathsToParameters (filesin['charmm_par'])
            system              = ImportSystem       ( filesin['charmm_psf'] , isXPLOR = True, parameters = parameters )
            system.coordinates3 = ImportCoordinates3 ( filesin['coordinates'] )
            self.define_NBModel(_type = 1, system = system)        
        elif systype == 2:
            mmModel        = MMModelOPLS.WithParameterSet ( filesin['opls_folder'] )            
            system         = ImportSystem       ( filesin['coordinates'])
            system.DefineMMModel ( mmModel )
            self.define_NBModel(_type = 1, system = system)          
        elif systype == 3:
            system = ImportSystem (filesin['coordinates'])
            system.Summary()
            #print ('mmModel',system.mmModel)
            #print ('qcModel',system.qcModel)
            #print ('nbModel',system.nbModel)

        '''
        psystem['system']        =  system
        psystem['name']          =  name
        print('color_palette', self.color_palette_counter)
        psystem['color_palette'] =  COLOR_PALETTE[self.color_palette_counter]
        #'''

        if name: 
            self.append_system_to_pdynamo_session(system, name = name, working_folder = HOME)
        else   : 
            name = system.label
            self.append_system_to_pdynamo_session(system, name = name, working_folder = HOME)
        
        self.vm_session.main_session.update_gui_widgets()


    def append_system_to_pdynamo_session (self, system, name = None, working_folder = None):
        """ Function doc """
        
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
                  'color_palette' : None , # will be replaced by a dict
                  'fixed_table'   : []   ,
                  'vismol_objects': {}   ,
                  'selections'    : {}   ,
                  'working_folder': HOME , #is a default folder that will be called to store simulation results, trajectories and log files. It is changed which the user changes the folder to a new simulation 
                  'logfile_data'  : {}   ,
                  'step_counter'  : 0    , #is a process counter that will be added to the name of each process executed inside easyhybrid
                   }
        
        if name:
            pass
        else:
            name = system.label
        
        try:
            psystem['system_original_charges'] =  list(system.AtomicCharges())
        except:
            psystem['system_original_charges'] = []
        psystem['system']                  =  system
        psystem['name']                    =  name
        psystem['color_palette']           =  COLOR_PALETTE[self.color_palette_counter]
        psystem['id']                      =  self.counter
        self.systems[psystem['id']]        =  psystem 
        
        #print('color_palette', self.color_palette_counter)
        #self.systems_list.append(psystem)
        
        self.active_id   = self.counter  
        self.counter    += 1

        if self.color_palette_counter >= len(COLOR_PALETTE)-1:
            self.color_palette_counter = 0
        else:
            self.color_palette_counter += 1
            
        self.build_vismol_object_from_pDynamo_system (name = 'initial coordinates' )#psystem['system'].label)
        if self.vm_session.main_session.selection_list_window.visible:
            self.vm_session.main_session.refresh_system_liststore()
            self.vm_session.main_session.selection_list_window.update_window(system_names = True, coordinates = False,  selections = False)
    
    def get_bonds_from_pDynamo_system(self, safety = 0.5, system_id = False):
        
        if system_id:
            pass
        else:
            system_id = self.active_id
        
        self.systems[system_id]['system'].BondsFromCoordinates3(safety = safety)
        
        raw_bonds =self.systems[system_id]['system'].connectivity.bondIndices
        bonds = []
        for bond in raw_bonds:
            bonds.append(bond[0])
            bonds.append(bond[1])
        
        self.systems[system_id]['bonds'] = bonds #self.systems[self.active_id]['system'].connectivity.bondIndices
        return True
        
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

    def get_fixed_atoms_from_system (self, system):
        """ Function doc """
        if system.freeAtoms is None:
            pass
        
        else:
            freeAtoms = system.freeAtoms
            freeAtoms = Selection.FromIterable (freeAtoms)
            selection_fixed = freeAtoms.Complement( upperBound = len (system.atoms ) )
            return list(selection_fixed)
            
    def define_free_or_fixed_atoms_from_iterable (self, fixedlist = None):
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
            #self.systems[self.active_id]['selections']["fixed atoms"] = list(selection_fixed)
            
            
            self.add_a_new_item_to_selection_list (system_id = self.active_id, 
                                                    indexes = list(selection_fixed), 
                                                    name    = 'fixed atoms')

        self.refresh_qc_and_fixed_representations(QC_atoms = False)
        return True

    def  _check_ref_vismol_object_in_pdynamo_system (self, system_id = None):
        """ Function doc """
        if system_id:
            pass
        else:
            system_id = self.active_id
        
        
        if self.systems[system_id]['vismol_object']:
            pass
        else:
            keys = self.systems[system_id]['vismol_objects'].keys()
            keys = list(keys)
            key = keys[-1]
            self.systems[system_id]['vismol_object'] = self.systems[system_id]['vismol_objects'][key]
    
    def check_charge_fragmentation(self, correction = True):
        """ Function doc """

        mm_residue_table = {}
        qc_residue_table = self.systems[self.active_id]['qc_residue_table']
        system           = self.systems[self.active_id]['system']
        qc_charge        = 0.0
        
        '''Here we are going to arrange the atoms that are not in the QC part, 
        but are in the same residues as some atoms of the QC part.'''  

        self._check_ref_vismol_object_in_pdynamo_system()
        
        for res in self.systems[self.active_id]['vismol_object'].residues:
            
            if res.resi in qc_residue_table.keys():
                
                mm_residue_table[res.resi] = []
                
                for atom in res.atoms:
                    index_v =  atom.index-1
                    index_p =  system.atoms.items[index_v].index
                    index_p =  system.atoms.items[index_v].label
                    charge  =  system.mmState.charges[index_v]
                    resn    = res.resn 
                    atom.charge = system.mmState.charges[index_v]
                    
                    if index_v in qc_residue_table[res.resi]:
                        qc_charge += atom.charge
                        pass
                        #print (resn, res.resi, index_v, index_p, charge, True )
                    
                    else:
                        #print (resn, res.resi, index_v, index_p, charge, False)
                        mm_residue_table[res.resi].append(index_v)
                
                
                
                #print(atom.index, atom.atomicNumber, system.mmState.charges[idx],self.systems[self.active_id]['vismol_object'].atoms[idx].resn )
            
        #print('mm_residue_table',mm_residue_table)
        '''Here we are going to do a rescaling of the charges of atoms of 
        the MM part but that the residues do not add up to an entire charge.''' 
        
        
        
        for resi in mm_residue_table.keys():
            
            total = 0.0
            for index in mm_residue_table[resi]:
                pcharge = system.mmState.charges[index]
                total += pcharge
            
            rounded  = float(round(total))
            diff     = rounded - total
            size     = len(mm_residue_table[resi])
            
            if size > 0:
                fraction = diff/size
                print('residue: ', resi, 'charge fraction = ',fraction)
            
                for index in mm_residue_table[resi]:
                    system.mmState.charges[index] += fraction
                    #total += pcharge
            else:
                pass
        print('QC charge from selected atoms: ',round(qc_charge) )
        #for atom in self.systems[self.active_id]['vismol_object'].atoms:
        #    print( atom.index, atom.name, atom.charge)
        #print('Total charge: ', sum(system.mmState.charges))
    
    
    def remove_item_from_selection_list (self, system_id = None, indexes = [], name  = 'sele'):
        """ Function doc """
        if system_id:
            pass
        else:
            system_id = self.active_id
            
        self.systems[system_id]['selections'].pop(name)            
    
    def set_selection_name (self, system_id):
        """ Function doc """
        if 'sel_name_counter' in  self.systems[system_id]:
            pass
        else:
            self.systems[system_id]['sel_name_counter'] = 0
        
        name = 'sel_'+ str(self.systems[system_id]['sel_name_counter'])
        
        loop = True
        while loop:
            if name in self.systems[system_id]['selections']:
                name = 'sel_'+ str(self.systems[system_id]['sel_name_counter'])
                self.systems[system_id]['sel_name_counter'] += 1
            else:
                loop = False
                
        return name
  
            
            
    
    def add_a_new_item_to_selection_list (self, system_id = None, indexes = [], name  = None):
        """ Function doc """
        
        if 'selections' in self.systems[system_id].keys():
            pass
        else:
            self.systems[system_id]['selections'] = {}
        
        
        if name:
            pass
        else:
            name = self.set_selection_name(system_id)
            #print('name', name)
            #name = 'sel_'+ str(self.sel_name_counter)
            #self.sel_name_counter += 1
        
        
        if system_id:
            pass
        else:
            system_id = self.active_id
            
        

            
        
        self.systems[system_id]['selections'][name] = indexes    
        if self.vm_session.main_session.selection_list_window.visible:
            self.vm_session.main_session.selection_list_window.update_window()

    def define_a_new_QCModel (self, parameters = None):
        """ Function doc """
        
        #print(parameters)
        
        '''Here we have to rescue the original electrical charges of the 
        MM model. This is postulated because multiple associations of QC 
        regions can distort the charge distribution of some residues. '''
        
        #self.systems[self.active_id]['system_original_charges']
        #for charge in self.systems[self.active_id]['system_original_charges']
        
        electronicState = ElectronicState.WithOptions ( charge = parameters['charge'], multiplicity = parameters['multiplicity'])
        self.systems[self.active_id]['system'].electronicState = electronicState

        qcModel         = QCModelMNDO.WithOptions ( hamiltonian = parameters['method'])
        

        
        if self.systems[self.active_id]['qc_table'] :
            
            '''This function reschedules the remaining loads in the MM part. The 
            sum of the charges in the MM region must be an integer value!'''
            self.check_charge_fragmentation()
            '''----------------------------------------------------------------'''
            
            self.systems[self.active_id]['system'].DefineQCModel (qcModel, qcSelection = Selection.FromIterable ( self.systems[self.active_id]['qc_table']) )
            self.refresh_qc_and_fixed_representations()# static = False )
            
            #print('define NBModel = ', self.nbModel)
            self.systems[self.active_id]['system'].DefineNBModel ( NBModelCutOff.WithDefaults ( ) )
            #self.systems[self.active_id]['selections']["QC atoms"] = self.systems[self.active_id]['qc_table']
            
            self.add_a_new_item_to_selection_list (system_id = self.active_id, 
                                                     indexes = self.systems[self.active_id]['qc_table'], 
                                                     name    = 'QC atoms')
            
            #print(self.systems[self.active_id]['selections']["QC atoms"])
        else:
            self.systems[self.active_id]['system'].DefineQCModel (qcModel)
            self.refresh_qc_and_fixed_representations()
            self.add_a_new_item_to_selection_list (system_id = self.active_id, 
                                                     indexes = range(0,self.systems[self.active_id]['system'].atoms), 
                                                     name    = 'QC atoms')
                                                     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    def refresh_vobject_qc_and_fixed_representations (self         ,
                                                      visObj = None,    
                                                 fixed_atoms = True , 
                                                    QC_atoms = True , 
                                                      static = True ):
        """ Function doc """
        
        
        system_id = visObj.easyhybrid_system_id
        if fixed_atoms:
            #system_id = visObj.easyhybrid_system_id
            #print ("system_id", system_id)
            
            if self.systems[system_id]['system'].freeAtoms is None:
                #print('pass')
                pass
            else:
                if self.systems[system_id]['fixed_table'] == []:
                    freeAtoms = self.systems[system_id]['system'].freeAtoms
                    freeAtoms                             = Selection.FromIterable (freeAtoms)
                    selection_fixed                       = freeAtoms.Complement( upperBound = len (self.systems[system_id]['system'].atoms ) )
                    self.systems[self.active_id]['fixed_table'] = list(selection_fixed)
                            #self.systems[self.active_id]['system'].freeAtoms = selection_free
                        
                        
                #print('set_color_by_index')
                #try:
            indexes = np.array(self.systems[system_id]['fixed_table'], dtype=np.int32)    
            color   = np.array(self.fixed_color, dtype=np.float32)    
            
            self.vm_session.set_color_by_index(vismol_object = visObj , 
                                               indexes       = indexes, 
                                               color         = color)
        
                #except:
                #    pass
        
        
        
        if QC_atoms:
            if self.systems[system_id]['system'].qcModel:
                self.systems[system_id]['qc_table'] = list(self.systems[system_id]['system'].qcState.pureQCAtoms)
                
                self.vm_session.show_or_hide_by_object (_type = 'spheres',  vobject = visObj, selection_table = range(0, len(visObj.atoms)),  show = False )
                self.vm_session.show_or_hide_by_object (_type = 'spheres', vobject = visObj, selection_table = self.systems[system_id]['qc_table'] , show = True )

                if static:
                    #self.vm_session.show_or_hide_by_object (_type = 'sticks',  vobject = visObj, selection_table = range(0, len(visObj.atoms)),  show = False)
                    self.vm_session.change_attributes_for_selected_atoms( _type = 'sticks', atoms = visObj.atoms , show = False )
                    self.vm_session.show_or_hide_by_object (_type = 'sticks' , vobject = visObj, selection_table = self.systems[system_id]['qc_table'] , show = True )

                else:
                    pass
                    self.vm_session.show_or_hide_by_object (_type = 'dynamic_bonds' , vobject = visObj, selection_table = self.systems[system_id]['qc_table'] , show = True )




























    
    def refresh_qc_and_fixed_representations (self, _all = False, 
                                               system_id = None , 
                                                  visObj = None,    
                                             fixed_atoms = True , 
                                                QC_atoms = True , 
                                                  static = True ):
        """ 
                
        _all = True/False applies the "ball and stick" and "color fixed atoms" representation
         to all vobjects. Only being used in load - serialization file
        
        """
        
        if system_id:
            pass
        else:
            system_id = self.active_id

        '''
        This loop is assigning the color of the fixed atoms to all objects 
        belonging to the pdynamo project that is active. 
        '''
        
        if _all:
            
            for key, visObj in self.vm_session.vismol_objects_dic.items():
                
                self.refresh_vobject_qc_and_fixed_representations ( 
                                                                    visObj = visObj     ,    
                                                               fixed_atoms = fixed_atoms, 
                                                                  QC_atoms = QC_atoms   , 
                                                                    static = static     )
                
                #system_id = visObj.easyhybrid_system_id
                #
                #if self.systems[system_id]['system'].qcModel:
                #    #print('\n\n\n\n system_id', system_id, visObj.name, visObj.easyhybrid_system_id, visObj.active )
                #    # Here we have to hide all the sticks and spheres so that there is no confusion in the representation of the QC region
                #    self.vm_session.show_or_hide_by_object (_type = 'spheres',  vobject = visObj, selection_table = range(0, len(visObj.atoms)),  show = False )
                #    self.vm_session.show_or_hide_by_object (_type = 'spheres', vobject = visObj, selection_table = self.systems[system_id]['qc_table'] , show = True )
                #
                #    if static:
                #        #self.vm_session.show_or_hide_by_object (_type = 'sticks',  vobject = visObj, selection_table = range(0, len(visObj.atoms)),  show = False)
                #        self.vm_session.change_attributes_for_selected_atoms( _type = 'sticks', atoms = visObj.atoms , show = False )
                #        self.vm_session.show_or_hide_by_object (_type = 'sticks' , vobject = visObj, selection_table = self.systems[system_id]['qc_table'] , show = True )
                #
                #    else:
                #        pass
                #        self.vm_session.show_or_hide_by_object (_type = 'dynamic_bonds' , vobject = visObj, selection_table = self.systems[system_id]['qc_table'] , show = True )
                #else:
                #    pass
            

        else:
            #if self.systems[system_id]['system'].qcModel:

                #self.systems[system_id]['qc_table'] = list(self.systems[system_id]['system'].qcState.pureQCAtoms)               
                
            for key, visObj in self.vm_session.vismol_objects_dic.items():
                if visObj.easyhybrid_system_id == self.active_id:
                    self.refresh_vobject_qc_and_fixed_representations ( 
                                                                      visObj = visObj     ,    
                                                                 fixed_atoms = fixed_atoms, 
                                                                    QC_atoms = QC_atoms   , 
                                                                      static = static     )
                        
    def merge_systems (self, system1 = None, system2 =  None, label = 'Merged System', summary = True):
        """ Function doc """
        system  = MergeByAtom ( system1, system2 )
        system.label = label
        self.define_NBModel(_type = 1, system = system)
        #system.define_NBModel( self.nbModel )
        
        if summary:
            system.Summary ( )
        
        '''
        psystem = {
                  'id'            : 0      ,  # access number (same as the access key in the self.systems dictionary)
                  'name'          : label  ,
                  'system'        : system ,  # pdynamo system
                  
                  'vismol_object' : None ,  # Vismol object associated with the system -> is the object that will 
                                            # undergo changes when something new is requested by the interface, for example: show the QC region
                  'active'        : False, 
                  'bonds'         : None ,
                  'sequence'      : None ,
                  'qc_table'      : None ,
                  'color_palette' : None , # will be replaced by a dict
                  'fixed_table'   : []   ,
                   }
        '''
        
        self.append_system_to_pdynamo_session (system)
        
    def prune_system (self, selection = None, label = 'Pruned System', summary = True):
        """ Function doc """
        p_selection   = Selection.FromIterable ( selection )
        system        = PruneByAtom ( self.systems[self.active_id]['system'], p_selection )
        self.define_NBModel(_type = 1, system = system)
        system.label  = label        
        if summary:
            system.Summary ( )
            
        '''
        psystem = {
                  'id'            : 0      ,  # access number (same as the access key in the self.systems dictionary)
                  'name'          : label  ,
                  'system'        : system ,  # pdynamo system
                  
                  'vismol_object' : None ,  # Vismol object associated with the system -> is the object that will 
                                            # undergo changes when something new is requested by the interface, for example: show the QC region
                  'active'        : False, 
                  'bonds'         : None ,
                  'sequence'      : None ,
                  'qc_table'      : None ,
                  'color_palette' : None , # will be replaced by a dict
                  'fixed_table'   : []   ,
                   }
        '''
            
            
        self.append_system_to_pdynamo_session (system)

    def get_coordinates_from_vismol_object_to_pDynamo_system (self, vismol_object ):
        """ Function doc """
        
        #print('\n\n')
        print('Loading coordinates from', vismol_object.name)
        #print('\n\n')
        
        for i, atom in enumerate(vismol_object.atoms):
            xyz = atom.coords(frame = -1)
            self.systems[self.active_id]['system'].coordinates3[i][0] = xyz[0]
            self.systems[self.active_id]['system'].coordinates3[i][1] = xyz[1]
            self.systems[self.active_id]['system'].coordinates3[i][2] = xyz[2]
    
    def get_sequence_from_pDynamo_system (self, system_id = None):
        """ Function doc """
        
        if system_id:
            pass
        else:
            system_id = self.active_id
        
        self.systems[system_id]['sequence'] = getattr ( self.systems[system_id]['system'], "sequence", None )
        if  self.systems[system_id]['sequence'] is None: 
            self.systems[system_id]['sequence'] = Sequence.FromAtoms ( self.systems[system_id]['system'].atoms, 
                                                                                     componentLabel = "UNK.1" )
        return True

    def get_atom_coords_from_pdynamo_system (self, system = None,  atom = None):
        if system:
            xyz = system.coordinates3[atom.index]
        else:
            xyz = self.systems[self.active_id]['system'].coordinates3[atom.index]
        return [float(xyz[0]),float(xyz[1]), float(xyz[2])]

    def get_atom_info_from_pdynamo_atom_obj (self, system_id = None, atom = None):
        """
        It extracts the information from the atom object, 
        belonging to pdynamo, and organizes it in the form 
        of a list that will be delivered later to build the 
        vismolObj
        
        """

        if system_id:
            pass
        else:
            system_id = self.active_id

        entityLabel = atom.parent.parent.label
        useSegmentEntityLabels = False
        if useSegmentEntityLabels:
            chainID = ""
            segID   = entityLabel[0:4]
        else:
            chainID = entityLabel[0:1]
            segID   = ""

        

        resName, resSeq, iCode = self.systems[system_id]['sequence'].ParseLabel ( atom.parent.label, fields = 3 )
        ##print(atom.index, atom.label,resName, resSeq, iCode,chainID, segID,  atom.atomicNumber, atom.connections)#, xyz[0], xyz[1], xyz[2] )
        
        index        = atom.index
        at_name      = atom.label
        at_resi      = int(resSeq)
        at_resn      = resName
        at_ch        = chainID
        at_symbol    = ATOM_TYPES_BY_ATOMICNUMBER[atom.atomicNumber] # at.get_symbol(at_name)
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
 
    def build_vismol_object_from_pDynamo_system (self                       , 
                                                 name = 'a_new_vismol_obj'  ,
                                                 system_id            = None,
                                                 vismol_object_active = True,
                                                 autocenter = True          ,
                                                 refresh_qc_and_fixed = True,
                                                 
                                                 ):
        """ Function doc """
        #print('\n\n\ build_vismol_object_from_pDynamo_system 736:', system_id, name)
        
        if system_id is not None:
            pass
        else:
            system_id = self.active_id
       #print('\n\n\ build_vismol_object_from_pDynamo_system 753:', system_id, name)
        
        name = str(self.systems[system_id]['step_counter'])+' '+name
        self.systems[system_id]['step_counter'] += 1
        
        self.get_bonds_from_pDynamo_system(safety = self.pdynamo_distance_safety, system_id = system_id)
        self.get_sequence_from_pDynamo_system(system_id = system_id)
        frames = []

        atoms = []     
        frame = []
        
        for atom in self.systems[system_id]['system'].atoms.items:
            xyz = self.get_atom_coords_from_pdynamo_system (atom   = atom, system  = self.systems[system_id]['system'])
            frame.append(xyz[0])
            frame.append(xyz[1])
            frame.append(xyz[2])
            
            atoms.append(self.get_atom_info_from_pdynamo_atom_obj(atom   = atom, system_id = system_id))
        

        
        frame = np.array(frame, dtype=np.float32)
        name  = os.path.basename(name)
        
        vismol_object  = VismolObject.VismolObject(name                           = name                                          ,    
                                                   atoms                          = atoms                                         ,    
                                                   vm_session                     = self.vm_session                            ,    
                                                   bonds_pair_of_indexes          = self.systems[system_id]['bonds']         ,    
                                                   trajectory                     = [frame]                                       ,  
                                                   color_palette                  = self.systems[system_id]['color_palette'] ,
                                                   auto_find_bonded_and_nonbonded = False               )
        
        vismol_object.easyhybrid_system_id = self.systems[system_id]['id']
        vismol_object.set_model_matrix(self.vm_session.glwidget.vm_widget.model_mat)
        vismol_object.active = vismol_object_active
        vismol_object._get_center_of_mass(frame = 0)
        
        if self.systems[system_id]['system'].qcModel:
            sum_x = 0.0 
            sum_y = 0.0 
            sum_z = 0.0
            
            self.systems[system_id]['qc_table'] = list(self.systems[system_id]['system'].qcState.pureQCAtoms)
            total = 0
            
            for atom_index in self.systems[system_id]['qc_table']:
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
        
        #self._check_ref_vismol_object_in_pdynamo_system()
        self.systems[system_id]['vismol_object'] = vismol_object
        
        self.vm_session.add_vismol_object_to_vismol_session (pdynamo_session  = self,
                                                                #rep              = True, 
                                                                vismol_object    = vismol_object, 
                                                                autocenter       = autocenter)
        
        if refresh_qc_and_fixed:
            self.refresh_qc_and_fixed_representations(system_id = system_id, visObj = vismol_object) 

        self.vm_session.glwidget.vm_widget.center_on_coordinates(vismol_object, center)
        self.vm_session.main_session.update_gui_widgets()
        return vismol_object

    def selections (self, _centerAtom, _radius, _method):
        """ Function doc """
        
        
        #print (_centerAtom)
        #print (_radius)
        #print (_method)
        self._check_ref_vismol_object_in_pdynamo_system()
        vismol_object = self.systems[self.active_id]['vismol_object']
        
        atomref = AtomSelection.FromAtomPattern( self.systems[self.active_id]['system'], _centerAtom )
        core    = AtomSelection.Within(self.systems[self.active_id]['system'],
                                                                      atomref,
                                                                      _radius)
                                                                      
        #core    = AtomSelection.Complement(self.systems[self.active_id]['system'],core)                                                
                                                        
                                                                      
        
        #print( core )
        
        if _method ==0:
            core    = AtomSelection.ByComponent(self.systems[self.active_id]['system'],core)
            core    = list(core)
            self.vm_session.selections[self.vm_session.current_selection].selecting_by_indexes (vismol_object   = vismol_object, 
                                                                                                              indexes = core , 
                                                                                                              clear   = True )
        
        if _method == 1:
            core    = AtomSelection.ByComponent(self.systems[self.active_id]['system'],core)
            core    = list(core)
            #'''******************** invert ? **********************
            inverted = []
            for i in range(0, len(vismol_object.atoms)):
                if i in core:
                    pass
                else:
                    inverted.append(i)
            
            core =  inverted
            self.vm_session.selections[self.vm_session.current_selection].selecting_by_indexes (vismol_object = vismol_object, 
                                                                                                      indexes = core, 
                                                                                                      clear   = True )

        if _method == 2:
            self.vm_session.selections[self.vm_session.current_selection].selecting_by_indexes (vismol_object   = vismol_object, 
                                                                                                              indexes = core , 
                                                                                                              clear   = True )
    
    def charge_summary (self, system = None):
        """ Function doc """
        
        if system == None:
            system = self.systems[self.active_id]['system']
            #self.systems[self.active_id]['system']
        
        self._check_ref_vismol_object_in_pdynamo_system()
        for res in self.systems[self.active_id]['vismol_object'].residues:
            for atom in res.atoms:
                index_v =  atom.index-1
                index_p =  system.atoms.items[index_v].index
                index_p =  system.atoms.items[index_v].label
                charge  =  system.mmState.charges[index_v]
                resn    = res.resn 
                atom.charge = system.mmState.charges[index_v]
                
               #print (resn, res.resi, index_v, index_p, charge )
                #print(atom.index, atom.atomicNumber, system.mmState.charges[idx],self.systems[self.active_id]['vismol_object'].atoms[idx].resn )
            
        #for atom in self.systems[self.active_id]['vismol_object'].atoms:
           #print( atom.index, atom.name, atom.charge)
        print('Total charge: ', sum(system.mmState.charges))
            
            
        #for atom in system.atoms.items:
        #    idx = atom.index
        #    print(atom.index, atom.atomicNumber, system.mmState.charges[idx],self.systems[self.active_id]['vismol_object'].atoms[idx].resn )
      
    def get_energy (self):
        """ Function doc """
        self.systems[self.active_id]['system'].Summary( )
        energy = self.systems[self.active_id]['system'].Energy( )
        return energy


    def import_data (self, _type = 'pklfile', data = None, logfile = None, first = 0 , last = -1, stride = 1, system_id = 0, vobject = None, name = None):
        """ Function doc """
        if _type == 'pklfile':
            frame = ImportCoordinates3 ( data )
            frame = list(frame) 
            print(list(frame))

            if vobject:
                pass
            else:
                vobject = self.build_vismol_object_from_pDynamo_system (
                                                                   name                 = name  ,
                                                                   system_id            = system_id,
                                                                   vismol_object_active = True        ,
                                                                   autocenter           = True        ,
                                                                   refresh_qc_and_fixed = False)
                vobject.frames = []            
            
            frame = np.array(frame, dtype=np.float32)
            vobject.frames.append(frame)
            self.refresh_qc_and_fixed_representations(_all = False, 
                                                 system_id = system_id,
                                                 visObj    = vobject,
                                                 fixed_atoms = True,
                                                 QC_atoms    = True,
                                                 static      = True,
                                                 ) 
            
            
            
            
            
        
        
        elif _type == 'pklfolder':
            self.import_trajectory ( traj = data, 
                                    first = first , 
                                     last = last, 
                                   stride = stride, 
                                system_id = system_id, 
                                  vobject = vobject, 
                                     name = name)
            
        elif _type == 'pklfolder2D':
            #print('vobject',vobject )
            vobject = self.import_2D_trajectory (traj = data, logfile = logfile, system_id = system_id, vobject = vobject, name = name)
        
            if logfile:
                data = parse_2D_scan_logfile (logfile)
                base = os.path.basename(logfile)
                print('vobject',vobject )
                if 'logfile_data' in  self.systems[system_id].keys():
                    if vobject.index in self.systems[system_id]['logfile_data']:
                        self.systems[system_id]['logfile_data'][vobject.index].append([base, data])
                
                    else:
                        self.systems[system_id]['logfile_data'][vobject.index] = []
                        self.systems[system_id]['logfile_data'][vobject.index].append([base, data])
                else:
                    self.systems[system_id]['logfile_data'] = {}
                    self.systems[system_id]['logfile_data'][vobject.index] = []
                    self.systems[system_id]['logfile_data'][vobject.index].append([base, data])
            
                self.vm_session.main_session.PES_analysis_window.OpenWindow(vobject = vobject)
                #self.vm_session.main_session.PES_analysis_window.on_coordinates_combobox_change(None)
        
        elif _type == 'pdbfile':
            pass
        elif _type == 'pdbfolder':
            pass
        elif _type == 'dcd':
            pass
        elif _type == 'crd':
            pass
        elif _type == 'xyz':
            pass
        elif _type == 'mol2':
            pass
        elif _type == 'netcdf':
            pass
        
        elif _type == 'log_file':
            if logfile:
                data = parse_2D_scan_logfile (logfile)
                base = os.path.basename(logfile)

                if 'logfile_data' in  self.systems[system_id].keys():
                    if vobject.index in self.systems[system_id]['logfile_data']:
                        self.systems[system_id]['logfile_data'][vobject.index].append([base, data])
                
                    else:
                        self.systems[system_id]['logfile_data'][vobject.index] = []
                        self.systems[system_id]['logfile_data'][vobject.index].append([base, data])
                
                
                else:
                    self.systems[system_id]['logfile_data'] = {}
                    self.systems[system_id]['logfile_data'][vobject.index] = []
                    self.systems[system_id]['logfile_data'][vobject.index].append([base, data])
            self.vm_session.main_session.PES_analysis_window.OpenWindow(vobject = vobject)
            #self.vm_session.main_session.PES_analysis_window.on_coordinates_combobox_change(None)
            pass
        else:
            pass
             
             
             
             
             
             
             
    def parse_logfile (self, logfile):
        """ Function doc """
        
        data  =  open(logfile, 'r')
        lines =  data.readlines()
        xlist = []
        ylist = []
        zlist = []
        
        lastline = lines[-1].split()
        x_size = int(lastline[0])
        y_size = int(lastline[1])
        
        #zlist = [[None]*(x_size+1)]*(y_size+1)
        #
        #for line in lines[1:]:
        #    line2 = line.split()
        #    x = int(line2[0])
        #    y = int(line2[1])
        #    print(x,y, line2[-1])
        #    zlist[y][x] = float(line2[-1])        
        
        #print (zlist)
        
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors

        rows = y_size+1
        cols = x_size+1
         
        zlist = [[0]*cols for _ in range(rows)]
        print (zlist)
        
        for line in lines[1:]:
            line2 = line.split()
            x = int(line2[0])
            y = int(line2[1])
            #print(x,y, line2[-1])
            zlist[y][x] = float(line2[-1]) 
            
        self.Z = zlist 
        
        #self.Z = [[56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], 
        #          [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036]]
        self.Y  = range(0,len(self.Z))
        self.X  = range(0,len(self.Z[0]))
        
        
        #N = 10 
        #X, Y = np.mgrid[-3:3:complex(0, N), -2:2:complex(0, N)]
        #Z1 = np.exp(-X**2 - Y**2)
        #Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2)
        #Z = (Z1 - Z2) * 2
        #print(Z)
        
        #Y  = range(0,len(zlist))
        #X  = range(0,len(zlist[0]))
        #Z  = zlist 
        
        fig, ax = plt.subplots(1, 1)
        #pcm = ax[0].pcolormesh(X, Y, Z, norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=-1.0, vmax=1.0, ), cmap='RdBu_r', shading='nearest')
        #fig.colorbar(pcm, ax=ax[0], extend='both')

        pcm = ax.pcolormesh(self.X, self.Y, self.Z, cmap='jet', vmin=np.min(self.Z), shading='gouraud')
        fig.colorbar(pcm, ax=ax)#, extend='both')

        plt.show()
        
        
             
    def import_2D_trajectory (self, traj = None,logfile = None, system_id = 0, vobject = None, name = None):
        """ Function doc """
        frames = []
        frame  = []

        
        if traj:
            if vobject:
                pass
            else:
                vobject = self.build_vismol_object_from_pDynamo_system (
                                                                   name                 = name  ,
                                                                   system_id            = system_id,
                                                                   vismol_object_active = True        ,
                                                                   autocenter           = True        ,
                                                                   refresh_qc_and_fixed = False)
                vobject.frames = []
            
            files = os.listdir(traj)
            files = sorted(files)
            
            vobject.trajectory2D_xy_indexes = {}
            vobject.trajectory2D_f_indexes  = {}

            n = 0
            for _file in files:
                if _file[-3:] == 'pkl':
                    frame = ImportCoordinates3 ( os.path.join(traj, _file) )
                    frame = list(frame) 
                    #print(list(frame))
                    
                    x_y = _file[5:-4].split('_')
                    
                    vobject.trajectory2D_xy_indexes[(int(x_y[0]), int(x_y[1]))] = n
                    vobject.trajectory2D_f_indexes[n] = (int(x_y[0]), int(x_y[1]))
                    
                    
                    frame = np.array(frame, dtype=np.float32)
                    vobject.frames.append(frame)
                    n+=1
            
            
            self.refresh_qc_and_fixed_representations(_all = False, 
                                                     system_id = system_id,
                                                     visObj    = vobject,
                                                     fixed_atoms = True,
                                                     QC_atoms    = True,
                                                     static      = True,
                                                     ) 
            return vobject
        else:
            pass
        
        
        '''--------------------------------------------------------------------'''
        #if logfile:
        #    #data = parse_2D_scan_logfile (logfile = traj[:-5]+'log')
        #    data = parse_2D_scan_logfile (logfile)
        #    #base = os.path.basename(traj[:-5]+'log')
        #    base = os.path.basename(logfile)
        #    
        #    #print(self.systems[system_id])
        #    
        #    if 'logfile_data' in  self.systems[system_id].keys():
        #        if vobject.index in self.systems[system_id]['logfile_data']:
        #            self.systems[system_id]['logfile_data'][vobject.index].append([base, data])
        #    
        #        else:
        #            self.systems[system_id]['logfile_data'][vobject.index] = []
        #            self.systems[system_id]['logfile_data'][vobject.index].append([base, data])
        #    
        #    
        #    else:
        #        self.systems[system_id]['logfile_data'] = {}
        #        self.systems[system_id]['logfile_data'][vobject.index] = []
        #        self.systems[system_id]['logfile_data'][vobject.index].append([base, data])
        #        
        #    
        #    #data = parse_2D_scan_logfile ('/home/fernando/ScanTraj2.log')
        #    #base = os.path.basename('/home/fernando/ScanTraj2.log')
        #    #self.systems[system_id]['logfile_data'][vobject.index].append([base, data])
        '''--------------------------------------------------------------------'''

        


        
    
    def import_trajectory (self, traj = None, first = 0 , last = -1, stride = 1, system_id = 0, vobject = None, name = None):
        """ Function doc """
        
        frames = []
        frame  = []
        trajectory = ImportTrajectory ( traj, self.systems[system_id]['system'] )
        trajectory.ReadHeader ( )
        
        # . Loop over the frames in the trajectory.
        phi = []
        psi = []
        
        if vobject:
            pass
        else:
            vobject = self.build_vismol_object_from_pDynamo_system (
                                                               name                 = name  ,
                                                               system_id            = system_id,
                                                               vismol_object_active = True        ,
                                                               autocenter           = True        ,
                                                               refresh_qc_and_fixed = False)
            vobject.frames = []
        
        while trajectory.RestoreOwnerData ( ):
            frame = []
            for atom in self.systems[system_id]['system'].atoms.items:
                xyz = self.get_atom_coords_from_pdynamo_system (atom   = atom, system = self.systems[system_id]['system'])
                frame.append(xyz[0])
                frame.append(xyz[1])
                frame.append(xyz[2])
            frame = np.array(frame, dtype=np.float32)
            #frames.append(frame)
            vobject.frames.append(frame)

        # . Finish up.
        trajectory.ReadFooter ( )
        trajectory.Close ( )
        #return frames
        self.refresh_qc_and_fixed_representations(system_id = system_id)           


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
        
        self.build_vismol_object_from_pDynamo_system (name = 'geometry optimization', autocenter = True)

    #---------------------------------------------------------------------------------
    def run_simulation(self, _parametersList = None):
        '''
        bsname = base name of the folder where will be created the next
        '''
        _parametersList["active_system"] = self.systems[self.active_id]['system']
        run = Simulation(_parametersList)
        run.Execute()        
        self.build_vismol_object_from_pDynamo_system (name = 'new_geometry', autocenter = False)
        
#======================================================================================================
 





