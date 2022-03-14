#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  easyhybrid_pDynamo_selection.py
#  
#  Copyright 2022 Fernando <fernando@winter>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
#from GTKGUI.gtkWidgets.filechooser import FileChooser
#from easyhybrid.pDynamoMethods.pDynamo2Vismol import *
import gc
import os

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')




class ImportTrajectoryWindow:
    """ Class doc """
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/easyhybrid_import_trajectory_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('import_trajectory_window')
            self.window.set_title('Import Trajectory Window')
            self.window.set_keep_above(True)
            '''--------------------------------------------------------------------------------------------'''
            
            self.combobox_pdynamo_system = self.builder.get_object('combobox_pdynamo_system')
            renderer_text = Gtk.CellRendererText()
            self.combobox_pdynamo_system.add_attribute(renderer_text, "text", 0)

            self.system_liststore = Gtk.ListStore(str)
            
            names = [ ]
            for system_key in self.easyhybrid_main.pDynamo_session.systems.keys():
                sys = self.easyhybrid_main.pDynamo_session.systems[system_key]
                names.append(sys['name'] )
                print(sys['name'] )
            
            for name in names:
                self.system_liststore.append([name])
            
            self.combobox_pdynamo_system.set_model(self.system_liststore)
            self.combobox_pdynamo_system.set_active(0)
            
            #methods = [
            #    "Conjugate Gradient" ,
            #    "FIRE"               ,
            #    "L-BFGS"             ,
            #    "Steepest Descent"   ,
            #    ]
            #for method in methods:
            #    self.method_store.append([method])
            #    print (method)
            #
            #self.methods_combo = self.builder.get_object('combobox_geo_opt')
            #self.methods_combo.set_model(self.method_store)
            #self.methods_combo.connect("changed", self.on_name_combo_changed)
            #self.methods_combo.set_model(self.method_store)
            #
            #renderer_text = Gtk.CellRendererText()
            #self.methods_combo.pack_start(renderer_text, True)
            #self.methods_combo.add_attribute(renderer_text, "text", 0)
            #'''--------------------------------------------------------------------------------------------'''
            self.combox = self.builder.get_object('combobox_trajectory_type')
            self.combox.set_active(0)
            self.window.show_all()
            self.Visible  = True
    
    def import_data (self, button):
        """ Function doc """
        entry_name    = None
        idnum     = self.combobox_pdynamo_system.get_active()
        text      = self.combobox_pdynamo_system.get_active_text()
        
        print(idnum, text )
        #traj_log      = int  ( self.builder.get_object('entry_traj_log').get_text() )
        #logFrequency  = int  ( self.builder.get_object('entry_log_frequence').get_text())
        #max_int       = int  ( self.builder.get_object('entry_max_int').get_text()  )
        #rmsd_tol      = float( self.builder.get_object('entry_rmsd_tol').get_text() )
        #entry_name    =        self.builder.get_object('entry_name').get_text()  
        #
        #save_trajectory = self.builder.get_object('checkbox_save_traj').get_active() 
    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
    
    
    def __init__(self, main = None):
        """ Class initialiser """
        self.easyhybrid_main     = main
        self.Visible             =  False        


