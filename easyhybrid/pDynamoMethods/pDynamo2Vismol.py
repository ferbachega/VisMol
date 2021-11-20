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



def get_atom_coords_from_pdynamo_system (system, atom, frame = None):

    if frame:
        xyz = system.coordinates3[atom.index]
        #print(atom.index, atom.label, atom.atomicNumber, atom.connections, xyz[0], xyz[1], xyz[2] )
        frame.append(float(xyz[0]))
        frame.append(float(xyz[1]))
        frame.append(float(xyz[2]))
    return frame

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
    print(atom.index, atom.label,resName, resSeq, iCode,chainID, segID,  atom.atomicNumber, atom.connections)#, xyz[0], xyz[1], xyz[2] )
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
        #print(atom.index, atom.label, atom.atomicNumber, atom.connections, xyz[0], xyz[1], xyz[2] )
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
        print(atom.index, atom.label,resName, resSeq, iCode,chainID, segID,  atom.atomicNumber, atom.connections, xyz[0], xyz[1], xyz[2] )
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

class pDynamoSession:
    """ Class doc """
    
    def __init__ (self, vismolSession = None):
        """ Class initialiser """
        self.vismolSession           = vismolSession
                                     
        self.name                    = 'pDynamo project'
        self.system                  = None
        self.vismol_object           = None
                                     
        self.bonds                   = None
        self.sequence                = None
        self.pdynamo_distance_safety = 0.5
                                    
        self.selection_qc_table      = None # indexes
        self.selection_fixed_table   = []   # indexes
        self.vismol_selection_qc     = None # vismol obj 
        
        self.nbModel = NBModelCutOff.WithDefaults ( )
        self.fixed_color = [0.5, 0.5, 0.5]
    
        self.pdm_systems = {
                            0:None
                            }
    
    
    def load_a_new_pDynamo_system_from_dict (self, filesin, systype):
        """ Function doc """
        
        if systype == 0:
            self.system              = ImportSystem       ( filesin['amber_prmtop'] )
            self.system.coordinates3 = ImportCoordinates3 ( filesin['coordinates'] )
            self.define_NBModel(_type = 1)
            
        if systype == 1:
            parameters          = CHARMMParameterFileReader.PathsToParameters (filesin['charmm_par'])
            self.system              = ImportSystem       ( filesin['charmm_psf'] , isXPLOR = True, parameters = parameters )
            self.system.coordinates3 = ImportCoordinates3 ( filesin['coordinates'] )
            self.define_NBModel(_type = 1)
        
        if systype == 2:
            mmModel             = MMModelOPLS.WithParameterSet ( filesin['opls_folder'] )            
            self.system         = ImportSystem       ( filesin['coordinates'])
            self.system.DefineMMModel ( mmModel )
            self.define_NBModel(_type = 1)
            
        if systype == 3:
            self.system = ImportSystem (filesin['coordinates'])
            self.system.Summary()
            print ('mmModel',self.system.mmModel)
            print ('qcModel',self.system.qcModel)
            print ('nbModel',self.system.nbModel)
            self.nbModel = self.system.nbModel
            print ('self.nbModel',self.nbModel)
            
        self.get_bonds_from_pDynamo_system(safety = self.pdynamo_distance_safety)

    
    
    def load_a_new_pDynamo_system (self, filein = None):
        """ Function doc """
        self.system = ImportSystem (filein)
        self.system.Summary()
        
    def get_bonds_from_pDynamo_system(self, safety = 0.5):
        
        self.system.BondsFromCoordinates3(safety = safety)
        self.bonds = self.system.connectivity.bondIndices
        
    
    def build_vismol_object_from_pDynamo_system (self                       , 
                                                 name = 'a_new_vismol_obj'  ,
                                                 vismol_object_active = True,
                                                 autocenter = True          ,):
        """ Function doc """
        self.get_bonds_from_pDynamo_system(safety = self.pdynamo_distance_safety)
        self.get_sequence_from_pDynamo_system()
        frames          = []

        atoms = []     
        frame = []
        
        for atom in self.system.atoms.items:
            xyz = self.get_atom_coords_from_pdynamo_system (atom)
            frame.append(xyz[0])
            frame.append(xyz[1])
            frame.append(xyz[2])
            
            atoms.append(self.get_atom_info_from_pdynamo_atom_obj(atom))
        
        #print ('frame', frame)
        
        frame = np.array(frame, dtype=np.float32)
        name  = os.path.basename(name)
        
        vismol_object  = VismolObject.VismolObject(name                           = name                ,    
                                                   atoms                          = atoms               ,    
                                                   vismolSession                  = self.vismolSession ,    
                                                   bonds_pair_of_indexes          = self.bonds          ,    
                                                   trajectory                     = [frame]             ,    
                                                   auto_find_bonded_and_nonbonded = False               )
        
        vismol_object.set_model_matrix(self.vismolSession.glwidget.vm_widget.model_mat)
        vismol_object.active = vismol_object_active
        vismol_object._get_center_of_mass(frame = 0)
        
        if self.system.qcModel:
            sum_x = 0.0 
            sum_y = 0.0 
            sum_z = 0.0
            
            self.selection_qc_table_table = list(self.system.qcState.pureQCAtoms)
            total = 0
            
            for atom_index in self.selection_qc_table_table:
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
            

        
        
        
        
        self.vismolSession.add_vismol_object_to_vismol_session (rep = True, vismol_object = vismol_object, autocenter =  autocenter)
        
        #center_on_coordinates
        self.vismolSession.glwidget.vm_widget.center_on_coordinates(vismol_object, center)
        
        self.vismol_object = vismol_object
        self.refresh_qc_and_fixed_representations()        
        return vismol_object
   
    
    def get_sequence_from_pDynamo_system (self):
        """ Function doc """
        # . Get the sequence.
        self.sequence = getattr ( self.system, "sequence", None )
        if self.sequence is None: self.sequence = Sequence.FromAtoms ( self.system.atoms, 
                                                                       componentLabel = "UNK.1" )

    def get_atom_coords_from_pdynamo_system (self, atom = None):

        xyz = self.system.coordinates3[atom.index]
        return [float(xyz[0]),float(xyz[1]), float(xyz[2])]

    def get_atom_info_from_pdynamo_atom_obj (self, atom):
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

        

        resName, resSeq, iCode = self.sequence.ParseLabel ( atom.parent.label, fields = 3 )
        #print(atom.index, atom.label,resName, resSeq, iCode,chainID, segID,  atom.atomicNumber, atom.connections)#, xyz[0], xyz[1], xyz[2] )
        
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
    
    def define_NBModel (self, _type = 1 , parameters =  None):
        """ Function doc """
        
        if _type == 0:
            self.nbModel = NBModelFull.WithDefaults ( )
        
        elif _type == 1:
            self.nbModel = NBModelCutOff.WithDefaults ( )
        
        elif _type == 2:
            self.nbModel = NBModelORCA.WithDefaults ( )
        
        elif _type == 3:
            self.nbModel = NBModelDFTB.WithDefaults ( )
        
        
        self.system.DefineNBModel ( self.nbModel )
        self.system.Summary ( )
        

        
    def get_energy (self):
        """ Function doc """
        self.system.Summary( )
        energy = self.system.Energy( )
        
        return energy
        
        
    
    def define_free_or_fixed_atoms_from_iterable (self, fixedlist):
        """ Function doc """
        selection_fixed            = Selection.FromIterable (fixedlist)
        self.selection_fixed_table = list(selection_fixed)
        selection_free             = selection_fixed.Complement( upperBound = len ( self.system.atoms ) )
        
        self.system.freeAtoms = selection_free
        self.refresh_qc_and_fixed_representations()
        return True
            
    def define_a_new_qcmodel (self, parameters):
        """ Function doc """
        
        print(parameters)
        
        electronicState = ElectronicState.WithOptions ( charge = parameters['charge'], multiplicity = parameters['multiplicity'])
        self.system.electronicState = electronicState

        qcModel         = QCModelMNDO.WithOptions ( hamiltonian = parameters['method'])
        

        
        if self.selection_qc_table :
            self.system.DefineQCModel (qcModel, qcSelection = Selection.FromIterable ( self.selection_qc_table ) )
            self.refresh_qc_and_fixed_representations()
            
            print('define NBModel = ', self.nbModel)
            self.system.DefineNBModel ( self.nbModel )
            #self.define_NBModel (_type = 1)
        
        else:
            self.system.DefineQCModel (qcModel)
            self.refresh_qc_and_fixed_representations()
            #self.vismolSession.show_or_hide (_type = 'spheres', selection = None,  show = True )
            #self.vismolSession.show_or_hide (_type = 'sticks', selection  = None,  show = True )
        
        
        
        
    def refresh_qc_and_fixed_representations (self):
        """ Function doc >>> 
        list(molecule.qcState.boundaryAtoms) 
        list(molecule.qcState.pureQCAtoms)
        list(molecule.qcState.qcAtoms)
        """
        #if self.selection_fixed_table:
        print('\n\n\nselection_fixed_table', self.selection_fixed_table)
        self.vismolSession.set_color_by_index(vismol_object = self.vismol_object, 
                                              indexes       = self.selection_fixed_table, 
                                              color         = self.fixed_color)
        
        if self.system.qcModel:

            self.selection_qc_table = list(self.system.qcState.pureQCAtoms)
            #boundaryAtoms     = list(self.system.qcState.boundaryAtoms)
            
            vismol_object = self.vismol_object
            
            # Here we have to hide all the sticks and spheres so that there is no confusion in the representation of the QC region
            self.vismolSession.selections[self.vismolSession.current_selection].selecting_by_indexes (vismol_object = vismol_object, indexes = range(0, len(vismol_object.atoms)))
            self.vismolSession.show_or_hide_by_object (_type = 'spheres',  vobject = vismol_object, selection_table = range(0, len(vismol_object.atoms)),  show = False )
            self.vismolSession.show_or_hide_by_object (_type = 'sticks',   vobject = vismol_object, selection_table = range(0, len(vismol_object.atoms)),  show = False )
            
            self.vismolSession.selections[self.vismolSession.current_selection].unselecting_by_indexes (vismol_object = self.vismol_object, indexes = range(0, len(vismol_object.atoms)))
            
            self.vismolSession.selections[self.vismolSession.current_selection].selecting_by_indexes (vismol_object = vismol_object, indexes = self.selection_qc_table)
            #print('\n\n',self.vismolSession.selections[self.vismolSession.current_selection].selected_atoms)
            
            self.vismolSession.show_or_hide_by_object (_type = 'spheres', vobject = vismol_object, selection_table = self.selection_qc_table, show = True )
            self.vismolSession.show_or_hide_by_object (_type = 'sticks' , vobject = vismol_object, selection_table = self.selection_qc_table, show = True )
            self.vismolSession.selections[self.vismolSession.current_selection].unselecting_by_indexes (vismol_object = self.vismol_object, indexes = range(0, len(vismol_object.atoms)))
            #print('\n\n',self.vismolSession.selections[self.vismolSession.current_selection].selected_atoms)

        else:
            pass


    def run_ConjugateGradientMinimize_SystemGeometry (self):
        """ Function doc """
        ConjugateGradientMinimize_SystemGeometry ( self.system                 ,
                                                   logFrequency         =  1   ,
                                                   maximumIterations    = 40 ,
                                                   rmsGradientTolerance =  0.1 )
        
        self.build_vismol_object_from_pDynamo_system (name = 'a_new_vismol_obj', autocenter = False)
