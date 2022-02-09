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
		'''
		'''
		self.molecule 		=_system
		self.baseName       = _baseName
		self.keywords 		= _keyWords
		self.QCatoms  		= []
		self.QCcharge 		= 0
		self.multiplicity 	= 1
		self.gradVectors 	= []
		self.molinFile      = None 
		self.inputFile 		= None
		self.atomsDict		= {}

		self.keywords.append("GRAD QMMM")
		# ver se é assim que se pega as cargas ainda
		self.charges = self.molecule.mmState.charges
		
		for i in self.molecule.atoms.items:
            symbol = GetAtomicSymbol( i.atomicNumber )
            index  = i.index
            x      = self.system.coordinates3[i.index, 0]
            y      = self.system.coordinates3[i.index, 1]
            z      = self.system.coordinates3[i.index, 2]
            atoms_dic[index] = [ symbol,x,y,z,charges[index] ]


		self.QCatoms		= list(self.molecule.qcState.pureQCAtoms)
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
	def write_input(self,_chg,_mult):
		'''
		Write the input files and grad vectors file
		'''

		mol_file  = open("mol.in","w")
		mop_file  = open( self.baseName + ".mop", "w")
		molInText = "\n{} 0".format( len(self.QCatoms) )
		mop_text  = ""

		for _key in self.keywords:
			mop_text += _key
		mop_text+="\n\n\n"

		for i in self.QCatoms:
			mop_text += "{} {} 1 {} 1 {} 1\n".format(atoms_dic[i][0],atoms_dic[i][2],atoms_dic[i][3],atoms_dic[i][4])

		for i in self.QCatoms:
			molInText += "{} {} {} {} {}\n".format(atoms_dic[i][0],atoms_dic[i][2],atoms_dic[i][3],atoms_dic[i][4],self.gradVectors[i])

		mop_file.write(mop_text)
		mop_file.close()

		mol_file.write(molInText)
		mol_file.close()



#======================================================================================================#
#======================================END OF FILE=====================================================#
#======================================================================================================#