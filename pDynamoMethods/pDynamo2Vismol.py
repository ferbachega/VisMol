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


class pDynamo2Vismol:
    """ Class doc """
    
    def __init__ (self, vismol_session = None):
        """ Class initialiser """
        self.vismol_session = vismol_session
        
        molecule = ImportSystem (  "/home/fernando/programs/EasyHybrid3/Coords/xyz/cyclohexane_single_frame.xyz"  )
        molecule.DefineQCModel ( QCModelMNDO.WithOptions ( hamiltonian = "am1" ) )
        molecule.Summary ( )
        
        molecule.BondsFromCoordinates3()
        molecule.Summary ( )
        
        for atom in molecule.atoms.items:
            xyz = molecule.coordinates3[atom.index]
            print(atom.index, atom.label, atom.atomicNumber, atom.connections, xyz[0], xyz[1], xyz[2] )
            #print('{} {} {} {} {} {}'.format(atom.index, atom.label, atom.connections, xyz[0], xyz[1], xyz[2]))
        
        bonds = molecule.connectivity.bondIndices
        print (bonds)
        
pdynamo_obj = pDynamo2Vismol()
        
        
        
        
        
        
        
        
        
        
