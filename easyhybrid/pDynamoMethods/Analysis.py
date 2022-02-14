#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#FILE = Analysis.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#==============================================================================

import os, sys, glob
import numpy as np

import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import seaborn as sns
from matplotlib.ticker import MaxNLocator
import matplotlib.colors as colors
from matplotlib.colors import BoundaryNorm
from collections import Counter

from commonFunctions import *

from pBabel                    import *                                     
from pCore                     import *                                     
from pMolecule                 import *                   
from pScientific               import *                                                             
from pScientific.Statistics    import *
from pScientific.Arrays        import *
from pSimulation               import *

logs1D = "/home/igorchem/CCDIR/LactatoPaper/LDL/SCANs1D/logs"
logs2D = "/home/igorchem/CCDIR/LactatoPaper/LDL/"

#*********************************************************************
class EnergyAnalysis:
	'''
	'''
	#------------------------------------------------
	def __init__(self, x, y, _type="1D"):
		'''
		'''
		self.energies1D 	= []
		self.energiesMatrix = np.zeros( (y, x), dtype=float )
		self.multiple1Dplot = []
		self.RC1            = []
		self.RC2            = []
		self.dimensions     = 0
		self.nplots1D       = 0
		self.xlen 			= x
		self.ylen  			= y
		self.Type 			= _type
		self.labely         = ""
		self.baseName 		= ""
		self.identifiers    = []

		if self.ylen > 0:
			self.dimensions = 2
		else:
			self.dimensions = 1

	#================================================
	def ReadLog(self, _fileName):
		'''

		'''
		self.baseName = _fileName
		reading = open(_fileName,'r')

		i = 0 
		energyTmp = []
		#----------------------------------
		if 	self.Type == "1D":
			for line in reading:
				if i > 0:
					lns = line.split()
					self.RC1.append( float(lns[1] ) )
					energyTmp.append( float(lns[2]) )
					self.energies1D.append( float(lns[2]) )
				i += 1

			self.multiple1Dplot.append(energyTmp)
			self.labely = "Potential Energy (kJ/mol)"

		#-----------------------------------
		elif self.Type == "1DRef":
			pass

		#----------------------------------
		elif self.Type == "2D":
			for line in reading:
				if i > 0:
					lns = line.split()
					self.RC1.append( float(lns[2] ) )
					self.RC2.append( float(lns[3] ) )
					m = int( lns[0])				
					n = int( lns[1])				
					self.energiesMatrix[n][m] = float(lns[4]) 
				i += 1
		
		#----------------------------------
		elif self.Type == "2DRef":
			for line in reading:
				lns = line.split()
				m = int( lns[0])				
				n = int( lns[1])				
				self.energiesMatrix[n][m] = float(lns[2])

		#----------------------------------
		elif self.Type == "WHAM1D":
			pass
		#----------------------------------
		elif self.Type == "WHAM2D":
			pass
		#----------------------------------
		elif self.Type == "FE1D":
			pass
		#----------------------------------
		elif self.Type == "FE2D":
			pass
		#----------------------------------
		self.nplots1D += 1
	
	#================================================
	def ReadLogs(self,_folder):
		'''
		'''
		_path = os.path.join(_folder,"")
		logs = glob.glob( _path + "*.log" )
		for log in logs:
			self.ReadLog(log)
			self.identifiers.append( os.path.basename(log[:-4]) )

	#==============================================
	def NormalizeEnergies(self):
		'''
		'''
		#------------------------------------------
		if self.Type == "1D":
			Min = 0
			if self.nplots1D == 1:
				Min = self.energies1D[0]
				for i in range( len(self.energies1D) ):
					self.energies1D[i] = self.energies1D[i] - Min

			elif self.nplots1D > 2:
				for k in range( self.nplots1D ):
					Min = self.multiple1Dplot[k][0]
					for i in range(len(self.multiple1Dplot)):
						self.multiple1Dplot[k][i] = self.multiple1Dplot[k][i] - Min
		
		#------------------------------------------
		if self.Type == "2D":
			self.energiesMatrix = self.energiesMatrix - np.min(self.energiesMatrix.min)


	#===============================================
	def FES_HL_SMO(self, logAB, logSMO, logFE):
		'''
		Free energy surface from a combination of High level QC method PES and semiempirical free energy
		Parameters:
			logAB:
			logSMO:
			logFE:
		'''
		pass
	#===============================================
	def Plot1D(self, label,SHOW=False):
		'''
		'''
		self.NormalizeEnergies()
		plt.plot(self.RC1,self.energies1D,'.b-')
		plt.xlabel(label)
		plt.ylabel(self.labely)		
		#--------------------------------------------
		_pathOut = self.baseName + "_1DenergyPlot.png"
		i=0
		while os.path.exists( _pathOut ):
			i+=1
			_pathOut = self.baseName + "_1DenergyPlot#{}.png".format(i)	
		plt.savefig(self.baseName+"_1Denergy.png",dpi=1000)
		#---------------------------------------------
		if SHOW:
			plt.show()
		
	#===============================================
	def MultPlot1D(self,label,SHOW=False):
		'''
		'''
		#---------------------------------------------
		self.NormalizeEnergies()
		x = np.linspace(0, self.xlen, self.xlen)
		for i in range(self.nplots1D):
			plt.plot(x,self.multiple1Dplot[i],label=self.identifiers[i],linestyle="-.")
		#---------------------------------------------
		plt.xlabel(label)
		plt.ylabel(self.labely)
		plt.legend()
		plt.savefig(self.baseName+"_1DmultiplePlots.png",dpi=1000)
		#---------------------------------------------
		if SHOW:
			plt.show()		

	#===============================================
	def Plot2D(self,contourlines,crd1label,crd2label,_xlim=None,_ylim=None,SHOW=False):
		'''
		'''
		x = range(self.xlen)
		y = range(self.ylen)		
		#-----------------------------------------------------
		if len(self.RC1) > 0:
			print(self.RC1[0],self.RC1[-1],self.RC2[0],self.RC2[-1])
			X = np.linspace(self.RC1[0],self.RC1[-1],self.xlen)
			Y = np.linspace(self.RC2[0],self.RC2[-1],self.ylen)
		#------------------------------------------------------
		else:
			if _xlim == None:
				_xlim = [ 0, self.xlen ]
				_ylim = [ 0, self.ylen ]
			X = np.linspace(_xlim[0],_xlim[1],self.xlen)
			Y = np.linspace(_xlim[0],_xlim[1],self.ylen)
		#------------------------------------------------------
		z = self.energiesMatrix
		#------------------------------------------------------
		fig, (ax0) = plt.subplots(nrows=1)
		vmin=z.min()
		vmax=z.max()
		#------------------------------------------------------
		levels = MaxNLocator(nbins=100).tick_values( z.min(), z.max() )
		cmap = plt.get_cmap("jet")
		#------------------------------------------------------
		norm = BoundaryNorm(levels, ncolors=cmap.N,	clip=True)
		norm= colors.PowerNorm(gamma=1./2.)
		norm= colors.Normalize(vmin=vmin, vmax=vmax)
		#------------------------------------------------------
		im = ax0.pcolormesh(X,Y,z, cmap=cmap, norm=norm, shading = "gouraud")
		am = ax0.contour(X,Y,z,contourlines, colors='k')		
		ax0.clabel(am, inline=1, fontsize=8, fmt='%1.1f',colors="k")		
		cbar = fig.colorbar(im, ax=ax0)
		cbar.ax.tick_params()
		#---------------------------------------------
		# Set the tick labels font
		axis_font = {'fontname':'Michroma', 'size':14}
		for tick in (ax0.xaxis.get_major_ticks()):
			tick.label.set_fontname('Arial')
			tick.label.set_fontsize(14)

		for tick in (ax0.yaxis.get_major_ticks()):
			tick.label.set_fontname('Dejavu')
			tick.label.set_fontsize(14) 
		#---------------------------------------------				
		ax0.set_xlabel(crd1label, **axis_font)
		ax0.set_ylabel(crd2label, **axis_font)
		fig.tight_layout()
		plt.savefig(self.baseName+"_PES.png",dpi=1000)
		if SHOW:
			plt.show()

#=====================================================================
#*********************************************************************
class DistanceAnalysis:
	'''
	'''
	#-----------------------------------------
	def __init__(self,_trajFolder,_system,t_time):
		'''
		'''
		self.trajFolder = _trajFolder
		self.molecule   = _system
		self.RG         = []
		self.RMS        = []
		self.distances1 = []
		self.distances2 = []
		self.total_time = t_time

		self.rc1_MF     = 0.0
		self.rc2_MF     = 0.0
		self.rg_MF      = 0.0
		self.rms_MF     = 0.0
		self.energies   = []

		if self.trajFolder[-4:] == ".dcd":			
			Duplicate(self.trajFolder,self.trajFolder[:-4]+".ptGeo",self.molecule)
			self.trajFolder = self.trajFolder[:-4]+".ptGeo"

		self.trajectory = ImportTrajectory( self.trajFolder , self.molecule )
		self.trajectory.ReadHeader()
        	
	#=================================================
	def CalculateRG_RMSD(self):
		'''
		Get Radius of Gyration and Root Mean Square distance for the trajectory
		'''
		masses  = Array.FromIterable ( [ atom.mass for atom in self.molecule.atoms ] )
		crd3    = ImportCoordinates3( os.path.join(self.trajFolder,"frame0.pkl") )
		system  = AtomSelection.FromAtomPattern ( self.molecule, "*:*:CA" )
		#------------------------------------------------------------------------------
		# . Calculate the radius of gyration.
		rg0 = crd3.RadiusOfGyration(selection = system, weights = masses)
		#------------------------------------------------------------------------------
		# . Save the starting coordinates.
		reference3 = Clone(crd3)  
		#------------------------------------------------------------------------------
		n = []
		m = 0             
		#-------------------------------------------------------------------------------
		while self.trajectory.RestoreOwnerData ( ):
			self.molecule.coordinates3.Superimpose ( reference3, selection = system, weights = masses )
			self.RG.append  ( self.molecule.coordinates3.RadiusOfGyration( selection = system, weights = masses ) )
			self.RMS.append ( self.molecule.coordinates3.RootMeanSquareDeviation( reference3, selection = system, weights = masses ) )
			n.append(m)
			m+=1
		# . Set up the statistics calculations.        
		rgStatistics  = Statistics(self.RG)
		rmsStatistics = Statistics(self.RMS)
		#-------------------------------------------------------------------------------
		# . Save the results.        
		textLog = open( self.trajFolder+"_MDanalysis", "w" ) 
		#-------------------------------------------------------------------------------
		_Text = "rg0 rgMean rgSD rgMax rgMin\n"
		_Text += "{} {} {} {} {}\n".format(rg0,rgStatistics.mean,rgStatistics.standardDeviation,rgStatistics.maximum,rgStatistics.minimum )
		_Text += "rmsMean rmsSD rmsMax rmsMin\n"
		_Text += "{} {} {} {}\n".format(rmsStatistics.mean,rmsStatistics.standardDeviation,rmsStatistics.maximum,rmsStatistics.minimum )
		#-------------------------------------------------------------------------------
		_Text += "Frame RG RMS\n"
		for i in range(len(self.RG)):
			_Text += "{} {} {}\n".format(i,self.RG[i],self.RMS[i])
		#--------------------------------------------------------------------------------
       
	#=================================================
	def ExtractFrames(self):
		'''
		'''
		#extractframe based on most common RMS
		#extractframe base on most commno RG and RMS
		#extractframe on rc1 and rc2 

		self.rc1_MF = Counter(self.distances1).most_common(1)[0][0]
		self.rc2_MF = Counter(self.distances2).most_common(1)[0][0]
		self.rms_MF = Counter(self.RMS).most_common(1)[0][0]
		self.rg_MF = Counter(self.RG).most_common(1)[0][0]

		distold = abs(self.RMS[0] - self.rms_MF)
		disnew  = 0.0
		fn      = 0 
		for i in range( len(self.RMS) ):
			distnew = abs(self.RMS[i] - self.rms_MF)
			if distnew < distold:
				distold = distnew
				fn = i
		rmsSystem = copySystem(self.molecule)
		rmsSystem.coordinates3 = ImportCoordinates3( os.path.join(self.trajFolder,"frame{}.pkl".format(fn) ) )
		ExportSystem("mostFrequentRMS.pdb",rmsSystem)
		ExportSystem("mostFrequentRMS.pkl",rmsSystem)
		#=============================================
		if len(self.distances2) > 0:
			distoldRC1 = abs(self.distances1[0] - self.rc1_MF)
			distoldRC2 = abs(self.distances2[0] - self.rc2_MF)
			distold    = abs(distoldRC1-distoldRC2)
			distnew    = 0.0
			fn         = 0 
			for i in range( len(self.distances1) ):
				distnew = abs( abs(self.distances1[i] - self.rc1_MF) -  abs(self.distances2[i] - self.rc2_MF) )
				if distnew < distold:
					distold = distnew
					fn = i 
			dsSystem = copySystem(self.molecule)
			dsSystem.coordinates3 = ImportCoordinates3( os.path.join(self.trajFolder,"frame{}.pkl".format(fn) ) )
			ExportSystem("mostFrequentRC1RC2.pdb",rmsSystem)
			ExportSystem("mostFrequentRC1RC2.pkl",rmsSystem)
			print(self.distances1[fn], self.distances2[fn] )

	#=================================================
	def PlotRG_RMS(self,SHOW=False):
		'''
		Plot graphs for the variation of Radius of Gyration and RMSD
		'''
		fig1, (ax1) = plt.subplots(nrows=1)

		n = np.linspace( 0, self.total_time, len(self.RG) )
		plt.plot(n, self.RG)
		ax1.set_xlabel("Time (ps)")
		ax1.set_ylabel("Radius of Gyration $\AA$")
		plt.savefig( os.path.join( self.trajFolder,"analysis_mdRG.png") )
		if SHOW: 
			plt.show()

		fig2, (ax2) = plt.subplots(nrows=1)
		#--------------------------------------------------------------------------
		plt.plot(n, self.RMS)
		ax2.set_xlabel("Time (ps)")
		ax2.set_ylabel("RMSD $\AA$")
		plt.savefig( os.path.join( self.trajFolder,"analysis_mdRMSD.png") )
		if SHOW:
			plt.show()        
		#---------------------------------------------------------------------------
		g = sns.jointplot(x=self.RG,y=self.RMS,kind="kde",cmap="plasma",shade=True)
		g.set_axis_labels("Radius of Gyration $\AA$","RMSD $\AA$",fontsize=12)
		plt.savefig( os.path.join( self.trajFolder,"rg_rmsd_biplot.png") )
		if SHOW:
			plt.show()        
	#=========================================================================
	def DistancePlots(self,RCs,SHOW=False):
		'''
		Calculate distances for the indicated reaction coordinates.
		'''		
		if len(RCs) == 2:
 			while self.trajectory.RestoreOwnerData():
 				self.energies.append( self.molecule.Energy(log=None) )
 				if RCs[0].nAtoms == 3:
 					self.distances1.append( self.molecule.coordinates3.Distance(RCs[0].atoms[0], RCs[0].atoms[1]) - self.molecule.coordinates3.Distance(RCs[0].atoms[1], RCs[0].atoms[2]) )
 				elif RCs[0].nAtoms == 2:
 					self.distances1.append( self.molecule.coordinates3.Distance(RCs[0].atoms[0], RCs[0].atoms[1]) )
 				if RCs[1].nAtoms == 3:                    
 					self.distances2.append( self.molecule.coordinates3.Distance(RCs[1].atoms[0], RCs[1].atoms[1]) - self.molecule.coordinates3.Distance(RCs[1].atoms[1], RCs[1].atoms[2]) )
 				elif RCs[1].nAtoms == 2:
 					self.distances2.append( self.molecule.coordinates3.Distance(RCs[1].atoms[0], RCs[1].atoms[1]) )
 				
		if len(RCs) == 1:
			while self.trajectory.RestoreOwnerData():
				self.energies.append( self.molecule.Energy(log=None) )
				if RCs[0].nAtoms == 3:
					self.distances1.append( self.molecule.coordinates3.Distance( RCs[0].atoms[0], RCs[0].atoms[1]) - self.molecule.coordinates3.Distance(RCs[0].atoms[1], RCs[0].atoms[2]) )
				elif RCs[0].nAtoms == 2:
					self.distances1.append( self.molecule.coordinates3.Distance(RCs[0].atoms[0], RCs[0].atoms[1]))
		
		#------------------------------------------------------------------------
		# . Save the results.        
		textLog = open( self.trajFolder+"_MDdistAnalysis", "w" )         
		_Text = ""
		if len(RCs) > 1:
			_Text = "Frame RC1 RC2 Energy\n"
			for i in range( len(self.distances1) ):
				_Text += "{} {} {} {}\n".format(i,self.distances1[i],self.distances2[i],self.energies[i])
		else:
			_Text = "Frame RC1 Energy\n"
			for i in range( len(rc1) ):
				_Text += "{} {} {}\n".format(i,self.distances1[i],self.energies[i])
		#-------------------------------------------------------------------------
		n = np.linspace( 0, self.total_time, len(self.distances1) )

		fig1, (ax1) = plt.subplots(nrows=1)
		plt.plot(n, self.energies)
		ax1.set_xlabel("Time (ps)")
		ax1.set_ylabel("Energy kJ/mol")
		plt.savefig(self.trajFolder+"_MDenergy.png")
		if SHOW:
			plt.show()
		#-------------------------------------------------------------------------
		fig2, (ax2) = plt.subplots(nrows=1)
		plt.plot(n, self.distances1,label=RCs[0].label,linestyle="-." )
		if len(RCs) ==2:
			plt.plot(n, self.distances2,label=RCs[1].label,linestyle="-.")
		plt.legend()
		ax2.set_xlabel("Time (ps)")
		ax2.set_ylabel("Distances $\AA$")			
		plt.savefig(self.trajFolder+"_MDdistAnalysis.png")
		if SHOW:
			plt.show()
		#-------------------------------------------------------------------------
		if len(RCs) == 2:
			g=sns.jointplot(x=self.distances1,y=self.distances2,kind="kde",cmap="plasma",shade=True)
			g.set_axis_labels(RCs[0].label,RCs[1].label)
			plt.savefig( os.path.join( self.trajFolder,"distanceBiplot.png") )
			if SHOW:
				plt.show()

#=================================================================================

