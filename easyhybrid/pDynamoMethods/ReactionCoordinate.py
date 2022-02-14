#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#FILE = ReactionCoordinate.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#==============================================================================

from commonFunctions import *
from pMolecule import *

#*****************************************************************************
class ReactionCoordinate:
	'''
	Class to set up and store reaction coordinate information
	'''
	def __init__(self,_atoms,_massConstraint,_type="Distance"):
		'''
		Types:
			distance
			multipleDistance
			Angle
			Dihedral
		'''
		self.atoms	        = _atoms
		self.nAtoms 		= len(_atoms)
		self.massConstraint = _massConstraint
		self.Type 			= _type
		self.weight13 		=  1.0
		self.weight31 		= -1.0
		self.period 		= 360.0
		self.increment      = 0.0
		self.minimumD  		= 0.0
		self.label 			= "Reaction Coordinate"

		if self.Type == "Distance":
			if self.nAtoms == 3:
				self.Type == "multipleDistance"

	#==========================================================================================================
	def SetInformation(self,_molecule,_dincre):
		'''		
		'''	
		self.increment = _dincre
		self.molecule  = _molecule
		
		sequence = getattr( self.molecule, "sequence", None )
		
		if self.Type == "multipleDistance":
			A1 = self.molecule.atoms.items[ self.atoms[0] ]
			A2 = self.molecule.atoms.items[ self.atoms[1] ]
			A3 = self.molecule.atoms.items[ self.atoms[2] ]
			A1res = A1.parent.label.split(".")
			A2res = A2.parent.label.split(".")
			A3res = A3.parent.label.split(".")
			#resName1 
			self.label =  A1.label + "(" + A1res[0] + A1res[1] + ")-"
			self.label += A2.label + "(" + A2res[0] + A2res[1] + ")--"
			self.label += "--"
			self.label += A3.label + "(" + A2res[0] + A2res[1] + ") $\AA$"
            #.-------------------------------------------------
			if self.massConstraint:				
				#------------------------------------------------
				atomic_n1 = A1.atomicNumber
				atomic_n3 = A3.atomicNumber
				mass_a1 = GetAtomicMass(atomic_n1)
				mass_a3 = GetAtomicMass(atomic_n3)
				self.weight13 = mass_a1 /(mass_a1+mass_a3)
				self.weight31 = mass_a3 /(mass_a1+mass_a3)
				self.weight31 = self.sigma_a3_a1*-1
				dist_a1_a2 = self.molecule.coordinates3.Distance( self.atoms[0], self.atoms[1] )
				dist_a2_a3 = self.molecule.coordinates3.Distance( self.atoms[1], self.atoms[2] )
				self.minimumD = ( self.weight13 * dist_a1_a2) - ( self.weight31 * dist_a2_a3*-1)
				#------------------------------------------------				
            #.----------------------
			else:
				dist_a1_a2 = self.molecule.coordinates3.Distance( self.atoms[0], self.atoms[1] )
				dist_a2_a3 = self.molecule.coordinates3.Distance( self.atoms[1], self.atoms[2] )
				self.minimumD =  dist_a1_a2 - dist_a2_a3
        #.--------------------------       
		elif self.Type == "Distance":
			A1 = self.molecule.atoms.items[ self.atoms[0] ]
			A2 = self.molecule.atoms.items[ self.atoms[1] ]
			A1res = A1.parent.label.split(".")
			A2res = A2.parent.label.split(".")
			self.label =  A1.label + "(" + A1res[0] + A1res[1] + ")--"
			self.label += A2.label + "(" + A2res[0] + A2res[1] + ") $\AA$"			
			self.minimumD = self.molecule.coordinates3.Distance( self.atoms[0], self.atoms[1] )

#===================================================================================================================
