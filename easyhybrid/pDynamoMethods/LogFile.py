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
class LogFile:
    '''
    Class to create and handle Logfiles of pDynamo
    '''

    def __init__(self, _filePath ):
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

    #---------------------------------------------------------
    def inputLine(self, _lineText):
        '''
        Class method to insert lines in the text container.
        '''
        self.text += _lineText 
        self.text +="\n"

   #---------------------------------------------------------
    def separator(self):
        '''
        '''
        self.text += "===================================================\n"

    #---------------------------------------------------------
    def close(self):
        '''
        Class method to close the file object.
        '''
        self.end = timer()
        cputime = self.end - self.start
        print("Cpu time: " + str(cputime) )

        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
        self.separator()
        self.text += "Finishing at: " + dt_string + "\n"
        self.text += "Elapsed time: " + str(cputime) + "\n"
        self.separator()
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