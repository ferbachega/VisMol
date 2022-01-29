#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#FILE = RelaxedScan.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#==============================================================================


#*****************************************************************************
class ReactionCoordinate:
	'''
	Class to set up and store reaction coordinate information
	'''
	def __init__(self,_atoms,_massConstraint,_type):
		'''
		Types:
			distance
			multipleDistance
			Angle
			Dihedral
		'''
		self.atomsIndices 	= atoms
		self.nAtoms 		= len(atoms)
		self.massConstraint = _massConstraint
		self.Type 			= _type
		self.weight13 		=  1
		self.weight31 		= -1
		self.period 		= 360.0
		self.increment      = 0.0
		self.minimumD  		= 0.0

	#==========================================================================
	def SetInformation(self,_molecule,_dincre):
		'''
		
		'''	
		self.increment = _dincre

		if self.Type == "multipleDistance":
            #.----------------------
            if self.massConstraint:
                atomic_n1 = self.molecule.atoms.items[ self.atoms[ndim][0] ].atomicNumber
                atomic_n3 = self.molecule.atoms.items[ self.atoms[ndim][2] ].atomicNumber
                mass_a1 = GetAtomicMass(atomic_n1)
                mass_a3 = GetAtomicMass(atomic_n3)
                self.sigma_a1_a3[ndim] = mass_a1 /(mass_a1+mass_a3)
                self.sigma_a3_a1[ndim] = mass_a3 /(mass_a1+mass_a3)
                self.sigma_a3_a1[ndim] = self.sigma_a3_a1[ndim]*-1
                dist_a1_a2 = self.molecule.coordinates3.Distance( self.atoms[ndim][0], self.atoms[ndim][1] )
                dist_a2_a3 = self.molecule.coordinates3.Distance( self.atoms[ndim][1], self.atoms[ndim][2] )
                self.DMINIMUM[ndim] = ( self.sigma_a1_a3[ndim] * dist_a1_a2) - ( self.sigma_a3_a1[ndim] * dist_a2_a3*-1)

            #.----------------------
            else:
                dist_a1_a2 = self.molecule.coordinates3.Distance( self.atoms[ndim][0], self.atoms[ndim][1] )
                dist_a2_a3 = self.molecule.coordinates3.Distance( self.atoms[ndim][1], self.atoms[ndim][2] )
                self.DMINIMUM[ndim] =  dist_a1_a2 - dist_a2_a3
        #.--------------------------       
        elif self.Type == "Distance":
            self.DMINIMUM[ndim] = self.molecule.coordinates3.Distance( self.atoms[ndim][0], self.atoms[ndim][1] )

#===============================================================================
