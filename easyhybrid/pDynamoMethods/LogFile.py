#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = LofFile.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

from pCore import *
from datetime import datetime
from timeit import default_timer as timer
#*************************************************************

HEADER = '''
#-----------------------------------------------------------------------------#
#                                                                             #
#                                EasyHybrid 3.0                               #
#                   - A pDynamo3 graphical user interface -                   #
#                                                                             #
#-----------------------------------------------------------------------------#
#                                                                             #
#             visit: https://sites.google.com/site/EasyHybrid/                #
#                                                                             #
#                                                                             #
#   EasyHybrid team:                                                          #
#   - Fernando Bachega                                                        #
#   - Igor Barden                                                             #
#   - Luis Fernando S M Timmers                                               #
#   - Martin Field                                                            #
#   - Troy Wymore                                                             #
#                                                                             #
#   Cite this work as:                                                        #
#   J.F.R. Bachega, L.F.S.M. Timmers, L.Assirati, L.B. Bachega, M.J. Field,   #
#   T. Wymore. J. Comput. Chem. 2013, 34, 2190-2196. DOI: 10.1002/jcc.23346   #
#                                                                             #
#-----------------------------------------------------------------------------#
'''


class LogFile:
    '''
    Class to create and handle Logfiles of pDynamo
    '''
    #-----------------------------------------------------------------
    def __init__(self, _filePath):
        '''
        Class constructor.
        Opens the file and initialize the text variable.
        '''

        self.start = timer()
        self.end   = 0 

        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

        self.filePath   = _filePath
        self.text       = "Log File for Simulation project on pDynamo make by EasyHybrid3.0!\n"
        self.text       += "Starting at: " + dt_string + "\n"
        self.separator()

        self.fileObj    = open( self.filePath,"w")

    #===================================================================
    def inputLine(self, _lineText):
        '''
        Insert lines in the text container.
        '''
        self.text += _lineText 
        self.text +="\n"

    #======================================================================
    def separator(self):
        '''
        Include a separator in the log Text.
        '''
        self.text += "===================================================\n"

    #======================================================================
    def close(self):
        '''
        Write and close the file object.
        '''
        #--------------------------------------------------------
        self.end = timer()
        cputime = self.end - self.start
        print("Cpu time: " + str(cputime) )
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
        #--------------------------------------------------------
        self.separator()
        self.text += "Finishing at: " + dt_string + "\n"
        self.text += "Elapsed time: " + str(cputime) + "\n"
        self.separator()
        #--------------------------------------------------------
        self.fileObj.write(self.text)
        self.fileObj.close()
    
    #---------------------------------------------------------
    def get_log(self):
        '''
        Class object to return a TextLogFileWriter pDynamo instance to use in individual methods
        '''
        logObj = TextLogFileWriter.WithOptions(self.filePath)
        return(logObj)
    
#====================================================================================================
