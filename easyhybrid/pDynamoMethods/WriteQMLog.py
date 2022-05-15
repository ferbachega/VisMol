#LDL QC/MM simulations of its catalysis reaction path

#-------------------------------------------------------------------
import sys                                                                                
sys.path.append("/home/igorchem/VisMol/easyhybrid/pDynamoMethods") 
import pymp
#--------------------------------------------------------------------
import os, glob
from commonFunctions import *
from CoreInterface import *
import SimulationsPreset

#==================================================
class WriteQMLog:
	'''
	'''
	def __init__(self,_system,_outFile):
		'''
		'''
		self.scratch = _system.scratch
		self.outname = _outFile

		self.text    = "{}".format( self.scratch.energyTerms["Potential Energy"] )
		self.charges = _system.AtomicCharges()

		self.outFile = open(_outFile,"w")

	#==============================================
	def write(self):
		'''
		'''
		#fill energy terms
		#number of atoms
		#fill atom labels, cooordinartes and charges
		#number of orbitals
		#basis set information
		#obital energies and occupancies 
		norbitals   = self.scratch.orbitalsP.numberOrbitals
		occupancies = self.scratch.orbitalsP.occupancies
		energies    = self.scratch.orbitalsP.energies		
		#self.scratch.overlapEigenValues 
		#overlap matrix
		block       = self.scratch.oneElectronMatrix.block

		print(len(block))
		for i in range(len(block)):
			if i % 10 == 0:
				self.text += "\n"
			self.text += "{} ".format(block[i])

		print(block[0],block[1])
		#molecular orbitals 
		orbitals    = self.scratch.orbitalsP.orbitals
		print(orbitals[1][0])


		outFile = open(self.outname,"w")
		outFile.write(self.text)
		outFile.close()
	
	#==============================================
	


#==================================================
#Scans 1D with the broken nad and arg 106
#--------------------------------------------------
#Setting reaction coordinates


