#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#FILE = ReactionCoordinate.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#==============================================================================

from commonFunctions import *

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

		if self.Type == "Distance":
			if self.nAtoms == 3:
				self.Type == "multipleDistance"

	#==========================================================================================================
	def SetInformation(self,_molecule,_dincre):
		'''		
		'''	
		self.increment = _dincre
		self.molecule  = _molecule

		if self.Type == "multipleDistance":
            #.----------------------
			if self.massConstraint:
				atomic_n1 = self.molecule.atoms.items[ self.atoms[0] ].atomicNumber
				atomic_n3 = self.molecule.atoms.items[ self.atoms[2] ].atomicNumber
				mass_a1 = GetAtomicMass(atomic_n1)
				mass_a3 = GetAtomicMass(atomic_n3)
				self.weight13 = mass_a1 /(mass_a1+mass_a3)
				self.weight31 = mass_a3 /(mass_a1+mass_a3)
				self.weight31 = self.sigma_a3_a1*-1
				dist_a1_a2 = self.molecule.coordinates3.Distance( self.atoms[0], self.atoms[1] )
				dist_a2_a3 = self.molecule.coordinates3.Distance( self.atoms[1], self.atoms[2] )
				self.minimumD = ( self.weight13 * dist_a1_a2) - ( self.weight31 * dist_a2_a3*-1)

            #.----------------------
			else:
				dist_a1_a2 = self.molecule.coordinates3.Distance( self.atoms[0], self.atoms[1] )
				dist_a2_a3 = self.molecule.coordinates3.Distance( self.atoms[1], self.atoms[2] )
				self.minimumD =  dist_a1_a2 - dist_a2_a3
        #.--------------------------       
		elif self.Type == "Distance":
			self.minimumD = self.molecule.coordinates3.Distance( self.atoms[0], self.atoms[1] )

#===================================================================================================================
