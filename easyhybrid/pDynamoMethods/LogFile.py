#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = LofFile.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

from pCore import *


class LogFile:
    '''
    Class to create and handle Logfiles of pDynamo
    '''

    def __init__(self, _filePath ):
        '''
        Class constructor.
        Opens the file and initialize the text variable.
        Task : 
            Customize the first message
        '''
        self.filePath   = _filePath
        self.text       = "Log File for Simulation project on pDynamo make by EasyHybrid3.0!"
        self.fileObj    = open( self.filePath,"w")

    def input_line(self, _lineText):
        '''
        Class method to insert lines in the text container.
        '''
        self.text += _lineText 
        self.text +="\n"
    
    def write(self):
        '''
        Class method to write the stored text on the file object.
        '''
        self.fileObj.write(self.text)
    
    def close(self):
        '''
        Class method to close the file object.
        '''
        self.fileObj.close()

    def get_log(self):
        '''
        Class object to return a TextLogFileWriter pDynamo instance to use in individual methods
        '''
        logObj = TextLogFileWriter(self.filePath)
        return(logObj)

    def UnitTest(self):
        '''
        Class method to Write some unit tests if needed in future
        '''
        pass