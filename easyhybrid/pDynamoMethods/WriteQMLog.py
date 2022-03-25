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
	def __init__(self,_system):
		'''
		'''
		self.scratch = _system.scratch
		self.text    = "{}".format( self.scratch.energyTerms["Potential Energy"])
		self.charges = _system.AtomicCharges()


	#==============================================
	def write(self,out_file):
		'''
		'''
		#fill energy terms
		#number of atoms
		#fill atom labels, cooordinartes and charges
		#number of orbitals
		#basis set information
		#obital energies and occupancies 
		norbitals   = scratch.orbitalsP.numberOrbitals
		occupancies = scratch.orbitalsP.occupancies
		energies    = scratch.orbitalsP.energies
		#overlap matrix
		block       = scratch.oneElectronMatrix.block
		#molecular orbitals 
		orbitals    = scratch.orbitalsP.orbitals

#==================================================
#Scans 1D with the broken nad and arg 106
#--------------------------------------------------
#Setting reaction coordinates

LDL = SimulationProject("OMTest")
LDL.LoadSystemFromSavedProject("LDL.pkl")

LDL.cSystem.Energy()
scratch = LDL.cSystem.scratch
	

