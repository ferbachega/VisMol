#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = EnergyRefinement.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#================================
import pymp
from commonFunctions import *
from pMolecule import *
from pMolecule.QCModel import *
from MopacQCMMinput import *
import os, glob, sys, shutil
import numpy as np 
#================================
#**********************************************************
class EnergyRefinement:
	'''
	Energy calculations for a set of structures using several QC methods and or External softwares
	'''
	def __init__(self,_refSystem,_trajFolder,_outFolder,_dims,_chg,_mult):
		'''
		Default constructor.
		Parameters:
			_refSystem : reference molecular information; System pDynamo class instance
			_trajFolder: folder path of the structures; string or path
			_outFolder : folder path where the results will be written; string or path
			_dims      : reaction coordinates size; list of integers
			_chg       : reference QC region charge; integer
			_multi     : reference QC region multiplicity; integer
		'''
		self.molecule 	 = _refSystem
		self.trajFolder  = _trajFolder  
		self.pureQCAtoms = list(self.molecule.qcState.pureQCAtoms)
		self.RC1 	 	 = []
		self.RC2 		 = []
		self.rc1CoordName= []
		self.rc2CoordName= []
		self.restart     = False
		self.xlen        = _dims[0]
		self.ylen        = _dims[1]
		self.charge 	 = _chg
		self.multiplicity= _mult
		self.text 		 = ""
		self.methods 	 = []

		i = 0
		_basepath =  os.path.join( _outFolder + "EnergyRefinement")		
		self.baseName = _basepath
		if not os.path.exists(self.baseName):
			os.makedirs(self.baseName)
		_path = os.path.join( _trajFolder,"")
		self.fileLists   = glob.glob(_path + "*.pkl")
		if self.ylen  == 0:
			self.energiesArray = pymp.shared.array( (self.xlen) , dtype='float')
			self.indexArrayX   = pymp.shared.array( (self.xlen) , dtype='uint8')
			self.indexArrayY   = pymp.shared.array( (self.xlen) , dtype='uint8')
		else:
			self.energiesArray = pymp.shared.array( (self.ylen,self.xlen) , dtype='float')
			self.indexArrayX   = pymp.shared.array( (self.ylen,self.xlen) , dtype='uint8')
			self.indexArrayY   = pymp.shared.array( (self.ylen,self.xlen) , dtype='uint8')
		self.SMOenergies   = None
	
	#=====================================================================================
	def GetQCCharge(self):
		'''
		'''
		qc_charge=0.0
		mmCharges = self.molecule.energyModel.mmAtoms.AtomicCharges()
		for qcatom in self.pureQCAtoms:
			qc_charge += mmCharges[qcatom]
		return(qc_charge)
	#=====================================================================================
	def ChangeQCRegion(self,_centerAtom,_radius):
		'''
		Redefine QC selection from a given atomic coordinates with a certain radius
		Parameters:
			_centerAtom:
			_radius    :
		'''
        #-----------------------------------------------------------------------------
		atomref      	 = AtomSelection.FromAtomPattern( self.molecule, _centerAtom )
		newSelection  	 = AtomSelection.Within(self.molecule,atomref,_radius)
		newSelection 	 = AtomSelection.ByComponent(self.molecule,newSelection)
		newSystem    	 = PruneByAtom(self.molecule, Selection(newSelection) )
		self.charge  	 = self.GetQCCharge()
		self.pureQCAtoms = list(newSelection)
		qcModel          = self.molecule.qcModel
		#-------------------------------------------------------------------------------
		self.molecule.electronicState = ElectronicState.WithOptions( charge = self.charge, multiplicity = self.multiplicity )
		self.molecule.DefineQCModel(qcModel, qcSelection=Selection(self.pureQCAtoms) )
        #---------------------------------------------------
        
	#=====================================================================================
	def RunInternalSMO(self,_methods,_NmaxThreads):
		'''
		Run energy refinement with the semiempirical hamiltonians available wihthin pDynamo
		Parameters:
			_methods:     List of Hamiltoninas
			_NmaxThreads: Number of maximum threds to be used in the parallel section
		'''
		self.SMOenergies = {}
		self.methods 	 = _methods
		NBmodel 	 	 = self.molecule.nbModel

		for smo in _methods:
			if VerifyMNDOKey(smo):
				with pymp.Parallel(_NmaxThreads) as p:
					for i in p.range( len(self.fileLists) ):
						qcModel = QCModelMNDO.WithOptions( hamiltonian = smo )
						self.molecule.electronicState = ElectronicState.WithOptions( charge = self.charge, multiplicity = self.multiplicity )
						self.molecule.DefineQCModel( qcModel, qcSelection=Selection(self.pureQCAtoms) )
						self.molecule.DefineNBModel( NBmodel )
						self.molecule.coordinates3 = ImportCoordinates3( self.fileLists[i],log=None )
						lsFrames= GetFrameIndex(self.fileLists[i][:-4])						
						if self.ylen > 0:
							self.energiesArray[ lsFrames[1], lsFrames[0] ] = self.molecule.Energy(log=None)
							self.indexArrayX[ lsFrames[1], lsFrames[0] ] = lsFrames[0]
							self.indexArrayY[ lsFrames[1], lsFrames[0] ] = lsFrames[1]
						else:
							self.energiesArray[ lsFrames[0] ] = self.molecule.Energy(log=None)
							self.indexArrayX[ lsFrames[0] ] = lsFrames[0]				

				#-----------------------------------------
				if self.ylen > 0:
					self.SMOenergies[smo] = self.energiesArray
					self.energiesArray = pymp.shared.array( (self.ylen,self.xlen) , dtype='float')	
				else:
					self.SMOenergies[smo] = self.energiesArray
					self.energiesArray = pymp.shared.array( (self.xlen) , dtype='float')		
			else:
				continue		
	#====================================================
	def RunInternalDFT(self,_functional,_basis,_NmaxThreads):
		'''
		Run energy refinement with the semiempirical hamiltonians available wihthin pDynamo
		Parameters:
			_methods:     List of Hamiltoninas
			_NmaxThreads: Number of maximum threds to be used in the parallel section
		'''
		self.SMOenergies = {}		
		self.methods.append(_functional+"/"+_basis)

		converger      = DIISSCFConverger.WithOptions( densityTolerance = 1.0e-10, maximumIterations = 250 )
		gridIntegrator = DFTGridIntegrator.WithOptions( accuracy = DFTGridAccuracy.Medium, inCore = True )
		qcModel        = None

		if not _functional == "hf": qcModel = QCModelDFT.WithOptions(converger=converger,functional="hf", orbitalBasis=_basis)
		else                      : qcModel = QCModelDFT.WithOptions(converger=converger,functional=_functional,gridIntegrator=gridIntegrator, orbitalBasis=_basis,)
		
		self.molecule.electronicState = ElectronicState.WithOptions(charge = self.charge)
		self.molecule.DefineQCModel( qcModel, qcSelection=Selection(self.pureQCAtoms) )
		self.molecule.DefineNBModel( NBmodel )
		
		
		with pymp.Parallel(_NmaxThreads) as p:
			for i in p.range( len(self.fileLists) ):				
				self.molecule.coordinates3 = ImportCoordinates3( self.fileLists[i],log=None )
				lsFrames= GetFrameIndex(self.fileLists[i][:-4])						
				if self.ylen > 0:
					self.energiesArray[ lsFrames[1], lsFrames[0] ] = self.molecule.Energy(log=None)
					self.indexArrayX[ lsFrames[1], lsFrames[0] ] = lsFrames[0]
					self.indexArrayY[ lsFrames[1], lsFrames[0] ] = lsFrames[1]
				else:
					self.energiesArray[ lsFrames[0] ] = self.molecule.Energy(log=None)
					self.indexArrayX[ lsFrames[0] ] = lsFrames[0]	
			#-----------------------------------------
			if self.ylen > 0:
				self.SMOenergies[ self.methods[0] ] = self.energiesArray
				self.energiesArray = pymp.shared.array( (self.ylen,self.xlen), dtype='float')	
			else:
				self.SMOenergies[ self.methods[0] ] = self.energiesArray
				self.energiesArray = pymp.shared.array( (self.xlen), dtype='float')

	#====================================================
	def RunMopacSMO(self,_methods,_keyWords):
		'''
		Create input for Mopac with its available Hamiltonians enabling the QC(QM)/MM potential
		Parameters:
			_methods: List of hamiltonians available in MOPAC 
			_NmaxThreads: Number of maximum threds to be used in the parallel section
		'''
		self.SMOenergies = {}
		self.methods     = _methods
		NBmodel          = self.molecule.nbModel	
		_mopacKeys       = _keyWords
		
		for smo in _methods:
			for i in range( len(self.fileLists) ):
				mop = MopacQCMMinput(self.molecule,self.baseName,self.fileLists[i][:-4],_mopacKeys,smo)
				mop.CalculateGradVectors()
				mop.write_input(self.charge,self.multiplicity)
				mop.Execute()
				lsFrames= GetFrameIndex(self.fileLists[i][:-4])						
				if self.ylen > 0:
					self.energiesArray[ lsFrames[1], lsFrames[0] ] = mop.GetEnergy()
					self.indexArrayX[ lsFrames[1], lsFrames[0] ] = lsFrames[0]
					self.indexArrayY[ lsFrames[1], lsFrames[0] ] = lsFrames[1]
				else:
					self.energiesArray[ lsFrames[0] ] = mop.GetEnergy()
					self.indexArrayX[ lsFrames[0] ] = lsFrames[0]	
			#----------------
			if self.ylen > 0:
					self.SMOenergies[smo] = self.energiesArray
					self.energiesArray    = pymp.shared.array( (self.ylen,self.xlen) , dtype='float')	
			else:
				self.SMOenergies[smo] = self.energiesArray
				self.energiesArray    = pymp.shared.array( (self.xlen) , dtype='float')	


	#====================================================
	def RunDFTB(self,_NmaxThreads):
		'''
		Perform energy refinement using the interface available on the pDynamo with the DFTB+ software, enabling QC(QM)/MM potential.
		Parameters:
			_NmaxThreads: Number of maximum threds to be used in the parallel section
		'''
		self.methods.append("DFTB")

		for i in p.range(0, len(self.fileLists) ):				
			fle2 = os.path.basename(self.fileLists[i][:-4])
			_scratch = os.path.join(self.baseName, fle2)				
			if not os.path.exists(_scratch):
				os.makedirs(_scratch)
							
			self.molecule.electronicState = ElectronicState.WithOptions(charge       = self.charge 		, 
			                                                          	multiplicity = self.multiplicity )
			#-------------------------------------------------------------
			_QCmodel = QCModelDFTB.WithOptions( deleteJobFiles = False   ,
			                                	randomScratch  = True    ,
			                                 	scratch        = _scratch,
			                                 	skfPath        = skfPath ,
			                                 	useSCC         = True    )
			#-------------------------------------------------------------
			NBmodel = NBModelDFTB.WithDefaults()
			self.molecule.DefineQCModel( _QCmodel, qcSelection=Selection(self.pureQCAtoms) )
			self.cSystem.DefineNBModel( self.NBmodel ) # reseting the non-bonded model
			#--------------------------------------------------------------------
			self.molecule.qcModel.maximumSCCIterations=1200
			energy = self.cSystem.Energy()			
			self.molecule.coordinates3 = ImportCoordinates3( self.fileLists[i] )
				
			if self.ylen > 0:
				self.energiesArray[ lsFrames[1], lsFrames[0] ] = self.molecule.Energy()					
				self.indexArrayX[ lsFrames[1], lsFrames[0] ]   = lsFrames[1]
				self.indexArrayY[ lsFrames[1], lsFrames[0] ]   = lsFrames[0]
				tmpText = "{}".format(self.energiesArray[ lsFrames[1], lsFrames[0] ])
				tmpLog.write(tmpText)
				tmpLog.close()
			else:					
				self.energiesArray[ lsFrames[0] ] = self.molecule.Energy()
				self.indexArrayX[ lsFrames[0] ]   = lsFrames[0]
				tmpText = "{}".format(self.energiesArray[ lsFrames[1], lsFrames[0] ])
				tmpLog.write(tmpText)
				tmpLog.close()
		 
	#====================================================
	def SetRestart4Orca(self):
		'''
		Set the files to be run in the energy refinement for ORCA with the restart option.
		The function will read a files named frame*_**.eTmp written in the folder with the energy of the frame.
		If the Orca Refinement run terminate succesfully, these files will be removed and the entire log file will be written.
		'''
		_path = os.path.join(self.baseName,"") 
		tmpList = glob.glob(_path+"*.eTmp")
		for fle in tmpList:
			lf = GetFrameIndex(fle[:-5])
			File = open(fle,'r')
			energy = File.read()
			if self.xlen > 0:
				self.indexArrayX[lf[0],lf[1]] 	= lf[0]
				self.indexArrayY[lf[0],lf[1]] 	= lf[1]
				self.energiesArray[lf[0],lf[1]]	= float(energy)
			else:
				self.indexArrayX[lf[0]] 	= lf[0]
				self.indexArrayY[lf[1]] 	= lf[1]
				self.energiesArray[lf[0]]	= float(energy)

		#-------------------------
		#remove files from list that already were calculated
		for fle in reversed(self.fileLists):
			fle2 = os.path.basename(self.fileLists[:-4])
			_scratch = os.path.join(self.baseName, fle2)				
			if os.path.exists(_scratch):
				self.fileLists.remove(fle)

	#====================================================
	def RunORCA(self,_method,_base,_NmaxThreads,_restart=False):
		'''
		Perform energy refinement using the interface available on the pDynamo with the ORCA software, enabling QC(QM)/MM potential.
		Parameters:
			_NmaxThreads: Number of maximum threds to be used in the parallel section
		'''
		self.methods.append(_method+"/"+_base)
		self.restart = _restart		
		if self.restart:
			self.SetRestart4Orca()				
		#---------------------------------------------------------
		#Initiate parallel run
		with pymp.Parallel(_NmaxThreads) as p:
			#----------------------------------------
			#Initiate Loop			
			for i in p.range(0, len(self.fileLists) ):				
				fle2 = os.path.basename(self.fileLists[i][:-4])
				_scratch = os.path.join(self.baseName, fle2)				
				if not os.path.exists(_scratch):
					os.makedirs(_scratch)
				#----------------------------------------------
				lsFrames= GetFrameIndex(self.fileLists[i][:-4])
				#----------------------------------------------
				tmpPath = os.path.join( self.baseName,"frame{}_{}.eTmp".format(lsFrames[0],lsFrames[1]) )
				tmpLog  = open(tmpPath,'w')
				tmpText = ""
				#---------------------------------------------
				options =  "\n% output\n"
				options +=  "print [ p_mos ] 1\n"
				options +=  "print [ p_overlap ] 5\n"
				options +=  "end # output"
				#...............................................................................................
				self.molecule.electronicState = ElectronicState.WithOptions(charge       = self.charge 		, 
				                                                          	multiplicity = self.multiplicity )
				#...............................................................................................
				QCmodel = QCModelORCA.WithOptions( keywords        = [ _method, _base, options], 
				                                   deleteJobFiles  = False                     ,
				                                   scratch         =_scratch                   )
				#...............................................................................................
				NBmodel = NBModelORCA.WithDefaults()
				self.molecule.DefineQCModel( QCmodel , qcSelection=Selection(self.pureQCAtoms) )
				self.molecule.DefineNBModel( NBmodel)
				self.molecule.coordinates3 = ImportCoordinates3( self.fileLists[i] )
				#---------------------------------------------------------------------------
				if self.ylen > 0:
					self.energiesArray[ lsFrames[1], lsFrames[0] ] = self.molecule.Energy()					
					self.indexArrayX[ lsFrames[1], lsFrames[0] ]   = lsFrames[1]
					self.indexArrayY[ lsFrames[1], lsFrames[0] ]   = lsFrames[0]
					tmpText = "{}".format(self.energiesArray[ lsFrames[1], lsFrames[0] ])
					tmpLog.write(tmpText)
					tmpLog.close()
				else:					
					self.energiesArray[ lsFrames[0] ] = self.molecule.Energy()
					self.indexArrayX[ lsFrames[0] ]   = lsFrames[0]
					tmpText = "{}".format(self.energiesArray[ lsFrames[1], lsFrames[0] ])
					tmpLog.write(tmpText)
					tmpLog.close()
		#--------------------
		self.TreatOrcaFiles()		
	#====================================================
	def TreatOrcaFiles(self):
		'''
		Rename orca files on the scratch folder, bringing them to the base folder with the name related with the respective frames
		'''
		_path    = os.path.join( self.baseName, "" )
		outFiles = glob.glob( _path+"frame*" )
		for out in outFiles:
			tmpPath     = os.path.join( out, "orcaJob.log" )
			outBasename = os.path.basename( out )
			finalPath   = os.path.join( self.baseName, out + ".out"  )
			shutil.move(outBasename,finalBasename)

	#====================================================
	def WriteLog(self):
		'''
		Write calculate energies to file.
		'''
		if self.ylen > 0:
			for smo in self.methods:
				for i in range(self.xlen):
					for j in range(self.ylen):
						self.text +="{} {} {} {}\n".format(self.indexArrayX[ i, j ],self.indexArrayY[ i, j ], self.SMOenergies[smo][i,j] - self.SMOenergies[smo][0,0], smo)
		else:
			for smo in self.methods:
				for i in range(self.xlen):
					self.text +="{} {} {}\n".format(self.indexArrayX[i], self.SMOenergies[smo][i] - self.SMOenergies[smo][0], smo)

		#--------------------------------------------------------------
		_filename = os.path.join(self.baseName+".log")
		#----------------------------
		logFile = open(_filename,'w')
		logFile.write(self.text)
		logFile.close()
		#----------------------------
		filesTmp = glob.glob("*.eTmp")
		for ftpm in filesTmp:
			os.remove(ftpm)		

#==========================================================


