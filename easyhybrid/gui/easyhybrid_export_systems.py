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
import sys

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')

sys.path.append(os.path.join(VISMOL_HOME,"easyhybrid/gui"))

from geometry_optimization_window import FolderChooserButton


class ImportTrajectoryWindow:
    """ Class doc """
    def OpenWindow (self, sys_selected = None):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/easyhybrid_import_trajectory_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('import_trajectory_window')
            self.window.set_title('Import Trajectory Window')
            self.window.set_keep_above(True)
            '''--------------------------------------------------------------------------------------------'''
            
            #------------------------------------------------------------------------------------
            self.combobox_pdynamo_system = self.builder.get_object('combobox_pdynamo_system')
            renderer_text = Gtk.CellRendererText()
            self.combobox_pdynamo_system.add_attribute(renderer_text, "text", 0)

            self.system_liststore = Gtk.ListStore(str, int)
            
            names = [ ]
            for key , system in self.easyhybrid_main.pDynamo_session.systems.items():
                name = system['name']
                self.system_liststore.append([name, int(key)])

            self.combobox_pdynamo_system.set_model(self.system_liststore)
            
            if sys_selected:
                self.combobox_pdynamo_system.set_active(sys_selected)
            else:
                self.combobox_pdynamo_system.set_active(0)
            #------------------------------------------------------------------------------------
            
            
            
            
            
            #'''--------------------------------------------------------------------------------------------------
            self.folder_chooser_button = FolderChooserButton(main =  self.window)
            self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
            #'''--------------------------------------------------------------------------------------------------
            
            
            
            #------------------------------------------------------------------------------------
            self.combobox_starting_coordinates = self.builder.get_object('vobjects_combobox')
            self.combobox_starting_coordinates.set_model(self.starting_coords_liststore)
            
            renderer_text = Gtk.CellRendererText()
            self.combobox_starting_coordinates.pack_start(renderer_text, True)
            self.combobox_starting_coordinates.add_attribute(renderer_text, "text", 0)
            
            size = len(self.starting_coords_liststore)
            self.combobox_starting_coordinates.set_active(size-1)
            #------------------------------------------------------------------------------------
            
            
            #------------------------------------------------------------------------------------
            self.builder.get_object('vobjects_combobox').set_sensitive(False)
            #------------------------------------------------------------------------------------

            #'''--------------------------------------------------------------------------------------------'''
            self.combox = self.builder.get_object('combobox_trajectory_type')
            self.combox.connect("changed", self.on_name_combo_changed)

            self.combox.set_active(0)
            self.window.show_all()
            self.Visible  = True
    
    def combobox_pdynamo_system (self, widget):
        """ Function doc """
        tree_iter = self.combobox_pdynamo_system.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_pdynamo_system.get_model()
            name, sys_id = model[tree_iter][:2]
            print (name, sys_id)
            #name, vobject_id = model[tree_iter][:2]
        
        
        self.starting_coords_liststore = Gtk.ListStore(str, int)
        for key, vobject  in self.easyhybrid_main.vm_session.vismol_objects_dic.items():
            print(key, vobject)
            if vobject.easyhybrid_system_id == sys_id:
                self.starting_coords_liststore.append([vobject.name, key])
        
        
        self.combobox_starting_coordinates = self.builder.get_object('vobjects_combobox')
        self.combobox_starting_coordinates.set_model(self.starting_coords_liststore)

    
    def on_vobject_combo_changed (self, widget):
        '''this combobox has the reference to the starting coordinates of a simulation'''
        #combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
        tree_iter = self.combobox_pdynamo_system.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_pdynamo_system.get_model()
            name, vobject_id = model[tree_iter][:2]
            print (name, model[tree_iter][:2])
            #name, vobject_id = model[tree_iter][:2]
    
    
    def on_name_combo_changed (self, widget):
        """ Function doc """
        if  self.combox.get_active() == 0:
            self.folder_chooser_button.sel_type = 'folder'
        else:
            self.folder_chooser_button.sel_type = 'file'

    def on_radio_button_change (self, widget):
        """ Function doc """
        if self.builder.get_object('radiobutton_import_as_new_object').get_active():
            self.builder.get_object('vobjects_combobox').set_sensitive(False)
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(True)

        if self.builder.get_object('radiobutton_append_to_a_vobject').get_active():
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(False)
            self.builder.get_object('vobjects_combobox').set_sensitive(True)

    def import_data (self, button):
        """ Function doc """
        entry_name     = None
        idnum          = self.combobox_pdynamo_system.get_active()
        text           = self.combobox_pdynamo_system.get_active_text()
        forder_or_file = self.folder_chooser_button.get_folder()
        
        
        #-----------------------------------------------------------------------------
        tree_iter = self.combobox_pdynamo_system.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_pdynamo_system.get_model()
            _name, sys_selected = model[tree_iter][:2]
            print (_name, sys_selected)
        #-----------------------------------------------------------------------------
       
       
       

        
        
        print(idnum, text,  forder_or_file)
        
        
        
        if self.builder.get_object('radiobutton_import_as_new_object').get_active():
            name    = self.builder.get_object('entry_create_a_new_vobj').get_text()
            vobject = None
        
        else:
            name = None
            #-----------------------------------------------------------------------------
            tree_iter = self.combobox_starting_coordinates.get_active_iter()
            if tree_iter is not None:
                
                '''selecting the vismol object from the content that is in the combobox '''
                model = self.combobox_starting_coordinates.get_model()
                name, vobject_id = model[tree_iter][:2]
                vobject = self.easyhybrid_main.vm_session.vismol_objects_dic[vobject_id]
            #-----------------------------------------------------------------------------
        
        
        #traj = os.path.join ( '/home/fernando/', 'NewTrajectory.ptGeo')
        self.easyhybrid_main.pDynamo_session.import_trajectory ( traj         = forder_or_file, 
                                                                 #first        =  0, 
                                                                 #last         = -1, 
                                                                 #stride       =  1,
                                                                 sys_selected =  sys_selected, 
                                                                 vobject      = vobject, 
                                                                 name         = name
                                                                 )
        

    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
    
    
    def __init__(self, main = None):
        """ Class initialiser """
        self.easyhybrid_main     = main
        self.Visible             =  False        
        self.starting_coords_liststore = Gtk.ListStore(str, int)

