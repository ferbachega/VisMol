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
		self.multiple2Dplot = []
		self.RC1            = []
		self.RC2            = []
		self.dimensions     = 0
		self.nplots1D       = 0
		self.nplots2D 		= 0
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
		self.baseName = _fileName[:-4]
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
			oldMethod = "none"
			method    = "none"
			for line in reading:
				lns = line.split()
				if oldMethod == "none":
					oldMethod = lns[2]
				method = lns[2]
				if not method == oldMethod:
					self.multiple1Dplot.append(energyTmp)
					self.identifiers.append(oldMethod)
					oldMethod = method
					self.nplots1D += 1	
					energyTmp = []		
				energyTmp.append( float(lns[1]) )
				self.energies1D.append( float(lns[1]) )

			self.multiple1Dplot.append(energyTmp)
			self.identifiers.append(method)
			self.labely = "Potential Energy (kJ/mol)"

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
			oldMethod = "none"
			method    = "none"
			for line in reading:
				lns = line.split()
				if oldMethod == "none":
					oldMethod = lns[3]
				method = lns[3]
				if not method == oldMethod:
					self.multiple2Dplot.append(self.energiesMatrix)
					self.identifiers.append(oldMethod)
					oldMethod = method
					self.nplots2D += 1
					self.energiesMatrix = np.zeros( (self.ylen, self.xlen), dtype=float )
				m = int( lns[0])				
				n = int( lns[1])				
				self.energiesMatrix[n][m] = float(lns[2])
			
			self.multiple2Dplot.append(self.energiesMatrix)
			self.identifiers.append(method)
			self.nplots2D += 1
		#----------------------------------
		elif self.Type == "WHAM1D":
			MaX = 0.0
			for line in reading:
				lns = line.split()
				pmf = float(lns[1])
				if pmf > MaX:
					MaX = pmf
				self.RC1.append( float(lns[0]) )
				if lns[1] == "inf":
					self.energies1D.append( 43434.0000 )
				else:
					self.energies1D.append( float(lns[1]) )

			for i in range(len(self.energies1D)):
				if self.energies1D[i] == 43434.0000:
					self.energies1D[i] == MaX

		#----------------------------------
		elif self.Type == "WHAM2D":
			m = 0
			n = 0
			MaX = 0.0
			for line in reading:
				pmf = float(lns[2])
				if pmf > MaX:
					MaX = pmf
				lns = line.split()
				self.RC1.append( float(lns[0]) )
				self.RC2.append( float(lns[1]) )
				if lns[2] == "inf":
					self.energiesMatrix[m][n] = 43434.0000
				else:
					self.energiesMatrix[m][n] = pmf				
				i +=1
				n +=1 						
				if i % self.xlen == 0:
					m += 1
					n = 0
			
			for j in range(self.xlen):
				for i in range(self.ylen):
					if self.energiesMatrix[i][j] == 43434.0000:
						self.energiesMatrix[i][j] = MaX

		#----------------------------------
		elif self.Type == "FE1D":
			for line in reading:
				lns = line.split()
				if lns[1] == "inf":
					self.energies1D.append( 180.0 )
				else:
					self.energies1D.append( float(lns[1]) )				
		#----------------------------------
		elif self.Type == "FE2D":
			for line in reading:
				lns = line.split()
				m = int( lns[0])				
				n = int( lns[1])
				if lns[1] == "inf":
					self.energies1D.append( 180.0 )	
				else:		
					self.energiesMatrix[n][m] = float(lns[2]) 
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
		Normalize energy arrays
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
		if self.Type == "2D" or self.Type == "WHAM2D" or self.Type == "FE2D" or self.Type == "2DRef":
			self.energiesMatrix = self.energiesMatrix - np.min(self.energiesMatrix)
	#===============================================
	def FES_HL_SMO(self, logPES, logSMO, logFE):
		'''
		Free energy surface from a combination of High level QC method PES and semiempirical free energy
		Parameters:
			logPES:
			logSMO:
			logFE:
		'''
		pass
	#===============================================
	def Plot1D(self, label,SHOW=False):
		'''
		Plot one dimensional energy plot.
		'''
		self.NormalizeEnergies()
		_pathOut = self.baseName
		if self.Type == "FE1D":
			self.RC1 = np.linspace( 0,len(self.energies1D),len(self.energies1D) )
			self.labely = "Free Energy (kJ/mol)"
			_pathOut += "_1Dfree_energy.png"
		elif self.Type == "WHAM1D":
			#self.RC1 = np.linspace( 0,len(self.energies1D),len(self.energies1D) )
			self.labely = "Potential of Mean Field (kJ/mol)"
			_pathOut += "_1D_PMF.png"
		elif self.Type == "1D" and self.Type == "1DRef" :
			_pathOut += "_1Denergy.png"
		#--------------------------------------------
		plt.plot(self.RC1,self.energies1D,'.k-')
		plt.xlabel(label)
		plt.ylabel(self.labely)		
		#--------------------------------------------
		plt.savefig(self.baseName+"_1Denergy.png",dpi=1000)
		#---------------------------------------------
		if SHOW:
			plt.show()		
	#===============================================
	def MultPlot1D(self,label,SHOW=False):
		'''
		Plot one-dimensinal energy plot for several methods
		'''
		#---------------------------------------------
		self.NormalizeEnergies()
		x = np.linspace(0, self.xlen, self.xlen )
		
		for i in range(self.nplots1D):
			plt.plot(x,self.multiple1Dplot[i],label=self.identifiers[i])
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
		Plot contour plot for potential, free energy and potential of mean field
		'''			
		#-----------------------------------------------------
		self.NormalizeEnergies()
		if len(self.RC1) > 0:
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
		levels = MaxNLocator(nbins=20).tick_values( z.min(), z.max() )
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

		_method = ""
		if len(self.identifiers) > 0:
			_method = self.identifiers[-1]
		plotName = self.baseName + _method
		if self.Type == "WHAM2D":
			plotName += "_PMF2D.png" 
		elif self.Type == "FE2D":
			plotName += "_FE2D.png"
		elif self.Type == "2D":
			plotName += "_2DPES.png"
		elif self.Type == "2DRef":
			plotName == "_2DPES.png"
		plt.savefig(plotName,dpi=1000)
		if SHOW:
			plt.show()
	#----------------------------------------------------------------------------------------
	def MultPlot2D(self,contourlines,crd1label,crd2label,_xlim=None,_ylim=None,SHOW=False):
		'''
		'''
		for i in range(self.nplots2D):
			self.identifiers.append( self.identifiers[i] )
			self.energiesMatrix = self.multiple2Dplot[i]
			self.Plot2D(contourlines,crd1label,crd2label,_xlim=None,_ylim=None,SHOW=False)

#=====================================================================




