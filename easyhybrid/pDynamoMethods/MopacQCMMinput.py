#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = MopacQCMMinput.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

from commonFunctions import *

#************************************************************
class MopacQCMMinput:
	'''
	Class to set methods to creat inputs for run QC/MM in mopac
	'''
	#---------------------------------------------------------
	def __init__(self,_system,_baseName,_keyWords):
		self.molecule 		=_system
		self.baseName       = _baseName
		self.keywords 		= [] 
		self.QCatoms  		= []
		self.MMatoms  		= []
		self.BoundaryAtoms 	= []
		self.QCcharge 		= 0
		self.multiplicity 	= 1
		self.gradVectors 	= []
		self.molinFile      = None 
		self.inputFile 		= None
		self.atomsDict		= {}
		
		# ver se é assim que se pega as cargas ainda
		self.charges = self.molecule.mmState.charges
		
		for i in self.molecule.atoms.items:
            symbol = GetAtomicSymbol( i.atomicNumber )
            index  = i.index
            x      = self.system.coordinates3[i.index, 0]
            y      = self.system.coordinates3[i.index, 1]
            z      = self.system.coordinates3[i.index, 2]
            atoms_dic[index] = [symbol,x,y,z,charges [index]]


		self.QCatoms		= list(self.molecule.qcState.qcAtoms)
        self.BoundaryAtoms	= list(self.molecule.qcState.boundaryAtoms)

	#==================================================================
	def CalculateGradVectors(self):
		'''
		Calculate the grad vectors for the mol.in
		'''
		PHI = 0.0 
		distance = 0.0
		#----------------------------------
		for j in range(len(self.QCatoms)):
			for i in self.molecule.atoms.items:
				distance = self.molecule.coordinates3.Distance( i, self.QCatoms[j] )
				PHI 	+= self.charges[i]/ distance

			PHI *= 332
			self.gradVectors.append(PHI)
			PHI=0		
	
	#===================================================================
	def write_input(self):
		'''
		Write the input files and grad vectors file
		'''
		pass

#======================================================================================================#
#======================================END OF FILE=====================================================#
#======================================================================================================#