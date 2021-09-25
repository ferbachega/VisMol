import os
import time
import multiprocessing
import numpy as np
#import vModel.atom_types as at 
#import vModel.cDistances as cdist
from   vModel import VismolObject
from pprint import pprint

cpdef load_pdb_file (infile = None, gridsize = 3, vismolSession =  None, frames_only = False):
    """ Function doc 

    gridsize =

    The longest covalent bond I can find is the bismuth-iodine single bond.
    The order of bond lengths is single > double > triple.
    The largest atoms should form the longest covalent bonds. 
    So we look at atoms in the lower right corner of the Periodic Table.
    The most likely candidates are Pb, Bi, and I.
    The experimental bond lengths are:
    Bi-I = 281 pm; Pb-I = 279 pm; I-I = 266.5 pm.
    So the polar covalent Bi-I bond is the longest covalent measured so far.

    """
    #-------------------------------------------------------------------------------------------
    #                                P D B     P A R S E R 
    #-------------------------------------------------------------------------------------------
    at  =  vismolSession.vConfig.atom_types
    with open(infile, 'r') as pdb_file:
        pdbtext = pdb_file.read()
        if 'ENDMDL' in pdbtext:
            rawframes = pdbtext.split('ENDMDL')
        else:
            #pass
            rawframes = pdbtext.split('END')
            #print (rawframes)
        
        
        
        frames    = get_list_of_frames_from_pdb_rawframes (rawframes = rawframes)
        
        if frames_only:
            return frames
        else:
            pass
            
        n    = 0 
        atoms = []
        while atoms == []:
            atoms     = get_list_of_atoms_from_rawframe(rawframes[n], gridsize, at =  at)
            n += 1
        
        
    #-------------------------------------------------------------------------------------------
    
    name = os.path.basename(infile)
    vismol_object  = VismolObject.VismolObject(name        = name, 
                                               atoms       = atoms, 
                                               vismolSession   = vismolSession, 
                                               trajectory  = frames)
    '''
    #-------------------------------------------------------------------------------------------
    #                                Bonded and NB lists 
    #-------------------------------------------------------------------------------------------
    
    #atoms, bonds_full_indices, bonds_pair_of_indices, NB_indices_list = cdist.generete_full_NB_and_Bonded_lists(atoms)
    
    #-------------------------------------------------------------------------------------------
    #print (bonds_pair_of_indices, NB_indices_list )
    #for atom in atoms:
    #    pprint (atom[8])
    #-------------------------------------------------------------------------------------------
    #                         Building   V I S M O L    O B J
    #-------------------------------------------------------------------------------------------

    
    #vismol_object.non_bonded_atoms  = NB_indices_list
    #vismol_object._generate_atomtree_structure()
    #vismol_object._generate_atom_unique_color_id()
    #vismol_object.index_bonds       = bonds_full_indices
    #vismol_object.import_bonds(bonds_pair_of_indices)
	'''

    #-------------------------------------------------------------------------------------------
    return vismol_object
    
    
    
    
    
    
cpdef get_list_of_atoms_from_rawframe(rawframe, gridsize = 3, at =  None):
    """ Function doc 

ATOM      1  N   THR A   1      -1.820  24.919  -5.344  1.00  0.00           N  
ATOM      2  CA  THR A   1      -1.256  24.379  -4.074  1.00  0.00           C  
ATOM     61  CB  ILE A   4      -7.386  -0.466   0.343  1.00  0.00           C  
ATOM  16     HO2 GLC  1       2.188   -0.704  0.939   1.00  300.00          H 0.0000   
ATOM     27  NH1 ARG A   5      68.029  23.029  29.719  1.00 38.75           N1+

iCode = ""
( serial       ,
  name         ,
  altLoc       ,
  resName      ,
  chainID      ,
  resSeq       ,
  u00          ,
  u11          ,
  u22          ,
  u01          ,
  u02          ,
  u12          ,
  segID        ,
  atomicNumber ,
  formalCharge ) = self._ParseFixedFormatLine ( line                    ,
                                                (  6, 11, int  , None ) ,
                                                ( 12, 16, None , ""   ) ,
                                                ( 16, 17, None , ""   ) ,
                                                ( 17, 20, None , ""   ) ,
                                                ( 21, 22, None , ""   ) ,
                                                ( 22, 27, int  , None ) ,
                                                ( 28, 35, float, None ) ,
                                                ( 35, 42, float, None ) ,
                                                ( 42, 49, float, None ) ,
                                                ( 49, 56, float, None ) ,
                                                ( 56, 63, float, None ) ,
                                                ( 63, 70, float, None ) ,
                                                ( 72, 76, None , ""   ) ,
                                                ( 76, 78, PeriodicTable.AtomicNumberFromSymbol, -1 ) ,
                                                ( 78, 80, None , ""   ) )

    """

    pdb_file_lines  = rawframe.split('\n')   
    atoms           = []
    cdef int index           = 0
    for line in pdb_file_lines:
        if line[:4] == 'ATOM' or line[:6] == 'HETATM':
            #try:
            at_name    = line[12:16].strip()
            at_pos     = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            
            at_resi    = int(line[22:27])
            at_resn    = line[17:20].strip()
            at_ch      = line[21]             
            
            
            at_symbol  = line[70:]
            #print ('at_symbol raw ',at_symbol )
            at_symbol  = at_symbol.strip()
            #print ('at_symbol raw2 ',at_symbol )

            if at_symbol == 'MG':
                at_symbol = 'Mg'
            else:
                pass
            
            
            if at_symbol in at.ATOM_TYPES.keys():
                #print ('at_symbol if ',at_symbol )
                pass
            else:
                #print ('at_symbol else1 ',at_symbol )
                at_symbol  = at.get_symbol(at_name)
                #print ('at_symbol else2 ',at_symbol )

            #print('at_symbol',at_name, at_symbol) 
            
            at_occup   = float(line[54:60])   #occupancy
            at_bfactor = float(line[60:66])
            at_charge  = 0.0

            cov_rad  = at.get_cov_rad (at_symbol)
            gridpos  = [int(at_pos[0]/gridsize), int(at_pos[1]/gridsize), int(at_pos[2]/gridsize)]
            #ocupan   = float(line[54:60])
            #bfactor  = float(line[60:66])
            
                            #0      1        2        3       4        5        6       7       8       9       10          11        12      
            atoms.append([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])
            #atoms.append([index, at_name, cov_rad,  None, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])
            #print (index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge )
            index += 1
            #except:
            #    pass
                #print(line)
    #print('atoms:', atoms)
    return atoms





cpdef get_list_of_frames_from_pdb_rawframes (rawframes = None):
    """ Function doc """
    n_processor = multiprocessing.cpu_count()
    pool        = multiprocessing.Pool(n_processor)
    frames      = pool.map(get_pdb_frame_coordinates, rawframes)
    framesout   = [] 
    
    
    for frame in frames:
        if frame:
            frame = np.array(frame, dtype=np.float32)
            framesout.append(frame)
        else:
            pass

    return framesout





cpdef get_pdb_frame_coordinates (str frame):
    """ Function doc """
    #print ('\nstarting: parse_pdb - building atom list')
    #initial          = time.time()
    pdb_file_lines    = frame.split('\n')
    
    #cdef float *frame_coordinates
    frame_coordinates = []
    
    
    for line in pdb_file_lines:
        if line[:4] == 'ATOM' or line[:6] == 'HETATM':
            #at_name  = line[12:16].strip()
            frame_coordinates.append(float(line[30:38]))
            frame_coordinates.append(float(line[38:46]))
            frame_coordinates.append(float(line[46:54]))
            #at_pos   = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])

    #frame_coordinates = np.array(frame_coordinates, dtype=np.float32)

    if len(frame_coordinates) == 0:
        return False
    else:
        return frame_coordinates
