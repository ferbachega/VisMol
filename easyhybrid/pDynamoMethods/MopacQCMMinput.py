#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = MopacQCMMinput.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################



class MopacQCMMinput:
	'''
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
		self.atomsDict		{}
		
		self.charges = self.molecule.energyModel.mmAtoms.AtomicCharges()
		
		for i in self.molecule.atoms.items:
            symbol = PeriodicTable.Symbol( i.atomicNumber )
            index  = i.index
            x      = self.system.coordinates3[i.index, 0]
            y      = self.system.coordinates3[i.index, 1]
            z      = self.system.coordinates3[i.index, 2]
            atoms_dic[index] = [symbol,x,y,z,charges [index]]


		self.QCatoms		= list( self.molecule.energyModel.qcAtoms.QCAtomSelection() )      #
        self.BoundaryAtoms	= list( self.molecule.energyModel.qcAtoms.BoundaryAtomSelection() )

	#---------------------------------------------------------
	def CalculateGradVectors(self):
		'''

		'''
		
	
	#---------------------------------------------------------
	def write_input(self):
		'''
		'''


