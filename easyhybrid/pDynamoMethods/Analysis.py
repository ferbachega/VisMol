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


logs1D = "/home/igorchem/CCDIR/LactatoPaper/LDL/SCANs1D/logs"
logs2D = "/home/igorchem/CCDIR/LactatoPaper/LDL/"

#*********************************************************************
class EnergyAnalysis:
	'''
	'''
	#------------------------------------------------
	def __init__(self,x,y,_type="1D"):
		'''
		'''
		self.energies1D 	= []
		self.energiesMatrix = np.zeros( (y, x),dtype=float )
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
		elif self.Type == "WHAM1D":
			pass
		#----------------------------------
		elif self.Type == "WHAM2D":
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
	def Plot1D(self, label):
		'''
		'''
		self.NormalizeEnergies()
		plt.plot(self.RC1,self.energies1D,'.b-')
		plt.xlabel(label)
		plt.ylabel(self.labely)
		_pathOut = self.baseName + "_1DenergyPlot.png"
		i=0

		while os.path.exists( _pathOut ):
			i+=1
			_pathOut = self.baseName + "_1DenergyPlot#{}.png".format(i)

		plt.savefig(self.baseName+"_1Denergy.png",dpi=1000)
		plt.show()
		
	#===============================================
	def MultPlot1D(self,label):
		'''
		'''
		self.NormalizeEnergies()
		x = np.linspace(0, self.xlen, self.xlen)
		for i in range(self.nplots1D):
			plt.plot(x,self.multiple1Dplot[i],label=self.identifiers[i],linestyle="-.")

		plt.xlabel(label)
		plt.ylabel(self.labely)
		plt.legend()
		plt.savefig(self.baseName+"_1DmultiplePlots.png",dpi=1000)
		plt.show()
		

	#===============================================
	def plot2D(self,contourlines,crd1label,crd2label):
		'''
		'''
		x = range(self.xlen)
		y = range(self.ylen)
		
		X = np.linspace(self.RC1[0],self.RC1[-1],self.xlen)
		Y = np.linspace(self.RC2[0],self.RC2[-1],self.ylen)
		
		z = self.energiesMatrix

		fig, (ax0) = plt.subplots(nrows=1)
		#setting plot parameters
		vmin=z.min()
		vmax=z.max()
		levels = MaxNLocator(nbins=100).tick_values( z.min(), z.max() )
		cmap = plt.get_cmap("jet")

		norm = BoundaryNorm(levels, ncolors=cmap.N,	clip=True)

		norm= colors.PowerNorm(gamma=1./2.)
		norm= colors.Normalize(vmin=vmin, vmax=vmax)

		im = ax0.pcolormesh(X,Y,z, cmap=cmap, norm=norm, shading = "gouraud")
		am = ax0.contour(X,Y,z,contourlines, colors='k')
		
		ax0.clabel(am, inline=1, fontsize=8, fmt='%1.1f',colors="k")
		
		cbar = fig.colorbar(im, ax=ax0)
		cbar.ax.tick_params()

		FontSize = 14

		# Set the tick labels font
		axis_font = {'fontname':'Michroma', 'size':14}
		for tick in (ax0.xaxis.get_major_ticks()):
			tick.label.set_fontname('Arial')
			tick.label.set_fontsize(FontSize)

		for tick in (ax0.yaxis.get_major_ticks()):
			tick.label.set_fontname('Dejavu')
			tick.label.set_fontsize(FontSize) 

		
		
		ax0.set_xlabel(crd1label, **axis_font)
		#ax0.set_xlabel(rcoord2, **axis_font)
		ax0.set_ylabel(crd2label, **axis_font)


		fig.tight_layout()
		plt.savefig(self.baseName+"_PES.png",dpi=1000)
		plt.show()
		#am = ax0.contour(x,y,z,lspacing, colors=lcolor)

		
#*********************************************************************
class DistanceAnalysis:

	'''
	'''
	#-----------------------------------------------
	def __init__(self,_trajFolder,_system):
		'''
		'''
		pass
	#===============================================
	def ReadLog(self, _fileName,_type="1D"):
		'''
		'''
		pass
	#===============================================
	def Plot1D(self):
		'''
		'''
		pass
	#===============================================
	def plot2D(self):
		'''
		'''
		pass 
#*********************************************************************
if __name__ == "__main__":
	'''
	test = EnergyAnalysis(24,20,_type="2D")
	test.ReadLog( os.path.join(logs2D,"LDLScan2Drm1_sca2D.log" ) )
	test.plot2D(10,"C2(LAC)-H--C4N(NAD)","O(LAC)-H--NE2(HIS193)")
	'''


