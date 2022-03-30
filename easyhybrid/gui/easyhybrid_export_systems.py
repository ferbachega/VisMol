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


class ExportDataWindow:
    """ Class doc """
    def OpenWindow (self, sys_selected = None):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/export_system_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('export_data_window')
            self.window.set_title('Export Data Window')
            self.window.set_keep_above(True)
            '''--------------------------------------------------------------------------------------------'''
            
            
            '''--------------------------------------------------------------------------------------------'''
            self.format_store = Gtk.ListStore(str)
            self.formats = {
                       0 : 'pkl - pdynamo system',
                       1 : 'pkl - pdynamo coordinates',
                       2 : 'pdb' ,
                       3 : 'xyz' ,
                       4 : 'mol' ,
                       5 : 'mol2',
                       }
                       
            for key, _format in self.formats.items():
                self.format_store.append([_format])
                #print (format)
            self.combobox_fileformat = self.builder.get_object('combobox_fileformat')
            self.combobox_fileformat.set_model(self.format_store)
            #self.formats_combo.connect("changed", self.on_name_combo_changed)
            self.combobox_fileformat.set_model(self.format_store)
            
            renderer_text = Gtk.CellRendererText()
            self.combobox_fileformat.pack_start(renderer_text, True)
            self.combobox_fileformat.add_attribute(renderer_text, "text", 0)
            '''--------------------------------------------------------------------------------------------'''
            
            
            
            
            
            
            
            
            
            
            
            '''--------------------------------------------------------------------------------------------'''
            self.psystem_liststore = Gtk.ListStore(str,int)
            names = [ ]
            for key , system in self.easyhybrid_main.pDynamo_session.systems.items():
                name = system['name']
                print ([name, int(key)])
                self.psystem_liststore.append([name, int(key)])
            
            self.psystem_combo = self.builder.get_object('combobox_pdynamo_system')
            self.psystem_combo.set_model(self.psystem_liststore)
            #self.psystem_combo.connect("changed", self.on_name_combo_changed)
            self.psystem_combo.set_model(self.psystem_liststore)
            
            renderer_text = Gtk.CellRendererText()
            self.psystem_combo.pack_start(renderer_text, True)
            self.psystem_combo.add_attribute(renderer_text, "text", 0)
            if sys_selected:
                self.psystem_combo.set_active(sys_selected)
            else:
                self.psystem_combo.set_active(0)
            
            '''--------------------------------------------------------------------------------------------'''
            
            
                        

            
            
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
            self.folder_chooser_button.set_folder(HOME)
            #------------------------------------------------------------------------------------

            
            

            


            #self.builder.get_object('entry_first').set_sensitive(False)
            #self.builder.get_object('label_first').set_sensitive(False)
            #
            #self.builder.get_object('entry_stride').set_sensitive(False)
            #self.builder.get_object('label_stride').set_sensitive(False)

            self.combobox_fileformat.set_active(0)
            self.window.show_all()
            self.Visible  = True
    
    def combobox_pdynamo_system (self, widget):
        """ Function doc """
        #print('combobox_pdynamo_system aqui')
        tree_iter = self.psystem_combo.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.psystem_combo.get_model()
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
        self.combobox_starting_coordinates.set_active(0)
    
    def on_vobject_combo_changed (self, widget):
        '''this combobox has the reference to the starting coordinates of a simulation'''
        #combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
            print (name, model[tree_iter][:2])
            #name, vobject_id = model[tree_iter][:2]
        
        
        if len(self.easyhybrid_main.vm_session.vismol_objects_dic[vobject_id].frames) > 1:
            print(self.easyhybrid_main.vm_session.vismol_objects_dic[vobject_id].name,
                  len(self.easyhybrid_main.vm_session.vismol_objects_dic[vobject_id].frames),'True')
            self.is_single_frame = True
            
            if self.combobox_fileformat.get_active( ) == 0:
                self.builder.get_object('entry_first').set_sensitive(False)
                self.builder.get_object('label_first').set_sensitive(False)
                self.builder.get_object('entry_stride').set_sensitive(False)
                self.builder.get_object('label_stride').set_sensitive(False)
            else:
                self.builder.get_object('entry_first').set_sensitive(True)
                self.builder.get_object('label_first').set_sensitive(True)            
                self.builder.get_object('entry_stride').set_sensitive(True)
                self.builder.get_object('label_stride').set_sensitive(True)
            
            
        else:
            print(self.easyhybrid_main.vm_session.vismol_objects_dic[vobject_id].name,
                  len(self.easyhybrid_main.vm_session.vismol_objects_dic[vobject_id].frames),'False')
            
            self.is_single_frame = False
            self.builder.get_object('entry_first').set_sensitive(False)
            self.builder.get_object('label_first').set_sensitive(False)
            self.builder.get_object('entry_stride').set_sensitive(False)
            self.builder.get_object('label_stride').set_sensitive(False)
        
    
    def on_combobox_fileformat (self, widget):
        """ Function doc """
        print('on_combobox_fileformat')
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
            print (name, model[tree_iter][:2])
            
        if len(self.easyhybrid_main.vm_session.vismol_objects_dic[vobject_id].frames) > 1:
            print(self.easyhybrid_main.vm_session.vismol_objects_dic[vobject_id].name,
                  len(self.easyhybrid_main.vm_session.vismol_objects_dic[vobject_id].frames),'True')
            self.is_single_frame = True
            if self.combobox_fileformat.get_active( ) == 0:
                self.builder.get_object('entry_first').set_sensitive(False)
                self.builder.get_object('label_first').set_sensitive(False)
                self.builder.get_object('entry_stride').set_sensitive(False)
                self.builder.get_object('label_stride').set_sensitive(False)
            else:
                self.builder.get_object('entry_first').set_sensitive(True)
                self.builder.get_object('label_first').set_sensitive(True)            
                self.builder.get_object('entry_stride').set_sensitive(True)
                self.builder.get_object('label_stride').set_sensitive(True)
        else:
            print(self.easyhybrid_main.vm_session.vismol_objects_dic[vobject_id].name,
                  len(self.easyhybrid_main.vm_session.vismol_objects_dic[vobject_id].frames),'False')
            
            self.is_single_frame = False
            self.builder.get_object('entry_first').set_sensitive(False)
            self.builder.get_object('label_first').set_sensitive(False)
            self.builder.get_object('entry_stride').set_sensitive(False)
            self.builder.get_object('label_stride').set_sensitive(False)
   
   
   
   
    def on_name_combo_changed (self, widget):
        """ Function doc """
        if  self.combox.get_active() == 0:
            self.folder_chooser_button.sel_type = 'folder'
        else:
            self.folder_chooser_button.sel_type = 'file'

    def on_radio_button_change (self, widget):
        """ Function doc """
        if self.builder.get_object('radiobutton_singlefile').get_active():
            self.is_single_file = True
            #if self.is_single_frame:
            #    self.builder.get_object('entry_first').set_sensitive(False)
            #    self.builder.get_object('label_first').set_sensitive(False)
            #    
            #    self.builder.get_object('entry_stride').set_sensitive(False)
            #    self.builder.get_object('label_stride').set_sensitive(False)
            #else:
            #    self.builder.get_object('entry_first').set_sensitive(True)
            #    self.builder.get_object('label_first').set_sensitive(True)
            #    
            #    self.builder.get_object('entry_stride').set_sensitive(True)
            #    self.builder.get_object('label_stride').set_sensitive(True)
            
        else:
            self.is_single_file = False

            #self.builder.get_object('entry_first').set_sensitive(True)
            #self.builder.get_object('label_first').set_sensitive(True)
            #
            #self.builder.get_object('entry_stride').set_sensitive(True)
            #self.builder.get_object('label_stride').set_sensitive(True)
        
        #if self.builder.get_object('radiobutton_append_to_a_vobject').get_active():
        #    self.builder.get_object('entry_create_a_new_vobj').set_sensitive(False)
        #    self.builder.get_object('vobjects_combobox').set_sensitive(True)

    def export_data (self, button):
        """ Function doc """
        
        '''Getting the correct format'''
        _format  = self.combobox_fileformat.get_active()
        _format  = self.formats[_format]
        _format  = _format.split()
        _format  = _format[0]
        '''--------------------------------------------'''
        
        folder   = self.folder_chooser_button.get_folder()
        filename = self.builder.get_object('entry_filename').get_text()
        
        print (folder, filename, _format)
        
        
        tree_iter = self.psystem_combo.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.psystem_combo.get_model()
            name, sys_id = model[tree_iter][:2]
            print (name, sys_id)

        self.easyhybrid_main.pDynamo_session.export_system (sys_id, filename, folder, _format)

        
        
    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
    
    
    def __init__(self, main = None):
        """ Class initialiser """
        self.easyhybrid_main     = main
        self.Visible             =  False        
        self.starting_coords_liststore = Gtk.ListStore(str, int)
        self.is_single_frame  = True
        
