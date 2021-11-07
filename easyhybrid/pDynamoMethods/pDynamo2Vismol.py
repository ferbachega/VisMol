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
    
    def __init__ (self, vismol_session = None):
        """ Class initialiser """
        self.system         = None
        self.vismol_session = vismol_session
        
        self.bonds          = None
        self.sequence       = None

        #
        #molecule = ImportSystem (  "/home/fernando/programs/EasyHybrid3/Coords/xyz/cyclohexane_single_frame.xyz"  )
        #molecule.DefineQCModel ( QCModelMNDO.WithOptions ( hamiltonian = "am1" ) )
        #molecule.Summary ( )
        #
        #molecule.BondsFromCoordinates3()
        #molecule.Summary ( )
        #
        #for atom in molecule.atoms.items:
        #    xyz = molecule.coordinates3[atom.index]
        #    print(atom.index, atom.label, atom.atomicNumber, atom.connections, xyz[0], xyz[1], xyz[2] )
        #    #print('{} {} {} {} {} {}'.format(atom.index, atom.label, atom.connections, xyz[0], xyz[1], xyz[2]))
        #
        #bonds = molecule.connectivity.bondIndices
        #print (bonds)
    
    
    def load_a_new_pDynamo_system_from_dict (self, filesin, systype):
        """ Function doc """
        
        if systype == 0:
            self.system              = ImportSystem       ( filesin['amber_prmtop'] )
            self.system.coordinates3 = ImportCoordinates3 ( filesin['coordinates'] )
    
        if systype == 1:
            parameters          = CHARMMParameterFileReader.PathsToParameters (filesin['charmm_par'])
            self.system              = ImportSystem       ( filesin['charmm_psf'] , isXPLOR = True, parameters = parameters )
            self.system.coordinates3 = ImportCoordinates3 ( filesin['coordinates'] )
    
    
    
    
    def load_a_new_pDynamo_system (self, filein = None):
        """ Function doc """
        #print(filein)
        
        #mmModel = MMModelOPLS.WithParameterSet ( '/home/fernando/programs/pDynamo3/parameters/forceFields/opls/protein' )
        #nbModel = NBModelMonteCarlo.WithDefaults ( )

        self.system = ImportSystem (filein)
        
        #self.system.DefineMMModel ( mmModel )
        #self.system.DefineNBModel ( nbModel )
        
        self.system.Summary()
        #self.system.Energy()
        
    def get_bonds_from_pDynamo_system(self):
        
        self.system.BondsFromCoordinates3()
        self.bonds = self.system.connectivity.bondIndices
        
    
    def build_vismol_object_from_pDynamo_system (self, name = 'a_new_vismol_obj'):
        """ Function doc """
        self.get_bonds_from_pDynamo_system()
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
                                                   vismolSession                  = self.vismol_session ,    
                                                   bonds_pair_of_indexes          = self.bonds          ,    
                                                   trajectory                     = [frame]             ,    
                                                   auto_find_bonded_and_nonbonded = False               )
        
        
        return vismol_object
   
    
    def get_sequence_from_pDynamo_system (self):
        """ Function doc """
        # . Get the sequence.
        self.sequence = getattr ( self.system, "sequence", None )
        if self.sequence is None: self.sequence = Sequence.FromAtoms ( self.system.atoms, 
                                                                       componentLabel = "UNK.1" )

    def get_atom_coords_from_pdynamo_system (self, atom = None):

        xyz = self.system.coordinates3[atom.index]
        #print(atom.index, atom.label, atom.atomicNumber, atom.connections, xyz[0], xyz[1], xyz[2] )
        #frame.append(float(xyz[0]))
        #frame.append(float(xyz[1]))
        #frame.append(float(xyz[2]))
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
        at_symbol    = self.vismol_session.vConfig.ATOM_TYPES_BY_ATOMICNUMBER[atom.atomicNumber] # at.get_symbol(at_name)
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
    
    
    def define_a_new_qcmodel (self, hamiltonian = 'am1'):
        """ Function doc """
        self.system.DefineQCModel ( QCModelMNDO.WithOptions ( hamiltonian = hamiltonian ) )

    
        
#pdynamo_obj = pDynamoObject()
        
#pdynamo_obj.load_a_new_pDynamo_system(filein = '/home/fernando/programs/VisMol/examples/xyz/cyclohexane_single_frame.xyz')
        
        
        
        
        
        
        
        
