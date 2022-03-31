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


class PotentialEnergyScanWindow():
    
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/PES_scan_window.glade'))
            self.builder.connect_signals(self)
            #
            self.window = self.builder.get_object('pes_scan_window')
            self.window.set_title('PES Scan Window')
            self.window.set_keep_above(True)            
            
            self.box_reaction_coordinate2 =  self.builder.get_object('box_reaction_coordinate2')           
            #'''--------------------------------------------------------------------------------------------'''
            #'''
            self.method_store = Gtk.ListStore(str)
            
            methods = ["simple distance", "combined distance","angle", 'dihedral']
            
            for method in methods:
                self.method_store.append([method])
                print (method)
            
            self.combobox_reaction_coord1 = self.builder.get_object('combobox_reaction_coord1')
            self.combobox_reaction_coord1.set_model(self.method_store)
            #self.combobox_reaction_coord1.connect("changed", self.on_name_combo_changed)
            
            self.combobox_reaction_coord2 = self.builder.get_object('combobox_reaction_coord2')
            self.combobox_reaction_coord2.set_model(self.method_store)
            #
            renderer_text = Gtk.CellRendererText()
            self.combobox_reaction_coord1.pack_start(renderer_text, True)
            self.combobox_reaction_coord1.add_attribute(renderer_text, "text", 0)
            
            self.combobox_reaction_coord2.pack_start(renderer_text, True)
            self.combobox_reaction_coord2.add_attribute(renderer_text, "text", 0)
            
            self.combobox_reaction_coord1.set_active(0)
            self.combobox_reaction_coord2.set_active(0)
            #'''--------------------------------------------------------------------------------------------'''
            self.method_store = Gtk.ListStore(str)
            
            methods = [
                "Conjugate Gradient" ,
                "FIRE"               ,
                "L-BFGS"             ,
                "Steepest Descent"   ,
                ]
            
            for method in methods:
                self.method_store.append([method])
                print (method)
            
            self.methods_combo = self.builder.get_object('combobox_methods')
            self.methods_combo.set_model(self.method_store)
            #self.methods_combo.connect("changed", self.on_name_combo_changed)
            self.methods_combo.set_model(self.method_store)
            #
            renderer_text = Gtk.CellRendererText()
            self.methods_combo.pack_start(renderer_text, True)
            self.methods_combo.add_attribute(renderer_text, "text", 0)
            #'''--------------------------------------------------------------------------------------------'''
            self.methods_combo.set_active(0)

            # #'''--------------------------------------------------------------------------------------------'''
            # 
            # '''
            # self.temp_scale_option_store = Gtk.ListStore(str)
            # 
            # temp_scale_options = ["constant", "linear","exponential"]
            # 
            # for temp_scale_option in temp_scale_options:
            #     self.temp_scale_option_store.append([temp_scale_option])
            #     print (temp_scale_option)
            # 
            # self.temp_scale_options_combo = self.builder.get_object('temperature_scale_option_combobox')
            # self.temp_scale_options_combo.set_model(self.temp_scale_option_store)
            # #self.temp_scale_options_combo.connect("changed", self.on_name_combo_changed)
            # self.temp_scale_options_combo.set_model(self.temp_scale_option_store)
            # #
            # renderer_text = Gtk.CellRendererText()
            # self.temp_scale_options_combo.pack_start(renderer_text, True)
            # self.temp_scale_options_combo.add_attribute(renderer_text, "text", 0)
            # #'''--------------------------------------------------------------------------------------------'''
            # self.temp_scale_options_combo.set_active(0)
            # 
            # 
            # 
            # 
            # job_list_canvas = self.builder.get_object('job_list_canvas')
            # 
            # software_list = [
            #     ("heating"        , '2002', "300"),
            #     ("Equilibration"  , '2004', "300"),
            #     ("Data-collection", '2004', "300"),
            #     #("Netbeans", 1996, "Java"),
            #     #("Chrome", 2008, "C++"),
            #     #("Filezilla", 2001, "C++"),
            #     #("Bazaar", 2005, "Python"),
            #     #("Git", 2005, "C"),
            #     #("Linux Kernel", 1991, "C"),
            #     #("GCC", 1987, "C"),
            #     #("Frostwire", 2004, "Java"),
            # ]
            #             
            # # Creating the ListStore model
            # self.software_liststore = Gtk.ListStore(str, str, str)
            # for software_ref in software_list:
            #     self.software_liststore.append(list(software_ref))
            # self.current_filter_language = None
            # 
            # # creating the treeview, making it use the filter as a model, and adding the columns
            # self.treeview = Gtk.TreeView(model=self.software_liststore)
            # for i, column_title in enumerate(
            #     ["Job", "nSteps", "temp"]
            # ):
            #     renderer = Gtk.CellRendererText()
            #     column = Gtk.TreeViewColumn(column_title, renderer, text=i)
            #     self.treeview.append_column(column)
            # 
            # ## creating buttons to filter by programming language, and setting up their events
            # #self.buttons = list()
            # #for prog_language in ["Java", "C", "C++", "Python", "None"]:
            # #    button = Gtk.Button(label=prog_language)
            # #    self.buttons.append(button)
            # #    button.connect("clicked", self.on_selection_button_clicked)
            # 
            # # setting up the layout, putting the treeview in a scrollwindow, and the buttons in a row
            # self.scrollable_treelist = Gtk.ScrolledWindow()
            # self.scrollable_treelist.set_vexpand(True)
            # #self.grid.attach(self.scrollable_treelist, 0, 0, 8, 10)
            # #self.grid.attach_next_to(
            # #    self.buttons[0], self.scrollable_treelist, Gtk.PositionType.BOTTOM, 1, 1
            # #)
            # #for i, button in enumerate(self.buttons[1:]):
            # #    self.grid.attach_next_to(
            # #        button, self.buttons[i], Gtk.PositionType.RIGHT, 1, 1
            # #    )
            # self.scrollable_treelist.add(self.treeview)
            # 
            # job_list_canvas.add(self.scrollable_treelist)
            # '''
            # 
            # #
            # #
            # #
            # #'''--------------------------------------------------------------------------------------------'''
            # #self.combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
            # #self.starting_coords_liststore = Gtk.ListStore(str)
            # #starting_coords = []
            # ##self.easyhybrid_main.pDynamo_session
            # #
            # #
            # #for key, visObj in self.easyhybrid_main.vm_session.vismol_objects_dic.items():
            # #    print(visObj.name, visObj.easyhybrid_system_id, visObj.active)
            # #    if visObj.easyhybrid_system_id == self.easyhybrid_main.pDynamo_session.active_id:
            # #        starting_coords.append(visObj.name)
            # #
            # #for method in starting_coords:
            # #    self.starting_coords_liststore.append([method])
            # #    print (method)
            # ##self.starting_coords_combo = self.builder.get_object('combobox_geo_opt')
            # #self.combobox_starting_coordinates.set_model(self.starting_coords_liststore)
            # ##self.combobox_starting_coordinates.connect("changed", self.on_name_combo_changed)
            # #self.combobox_starting_coordinates.set_model(self.starting_coords_liststore)
            # #
            # #renderer_text = Gtk.CellRendererText()
            # #self.combobox_starting_coordinates.pack_start(renderer_text, True)
            # #self.combobox_starting_coordinates.add_attribute(renderer_text, "text", 0)
            # #
            # #size = len(starting_coords)
            # #self.combobox_starting_coordinates.set_active(size-1)
            # '''--------------------------------------------------------------------------------------------'''

            self.window.show_all()
            self.box_reaction_coordinate2.set_sensitive(False)
            self.Visible  = True
    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
    
    
    def __init__(self, main = None):
        """ Class initialiser """
        self.easyhybrid_main     = main
        self.Visible             =  False        
        self.residue_liststore = Gtk.ListStore(str, str, str)


    def add_job_to_list (self, button):
        """ Function doc """
        self.software_liststore.append(list(['aqui', '123' ,'cocozao']))

    def change_radio_button_reaction_coordinate (self, widget):
        """ Function doc """
        #radiobutton_bidimensional = self.builder.get_object('radiobutton_bidimensional')
        if self.builder.get_object('radiobutton_bidimensional').get_active():
            self.box_reaction_coordinate2.set_sensitive(True)
        else:
            self.box_reaction_coordinate2.set_sensitive(False)
        
        #print(widget)

    #-------------------------------------------------------------------------------
    def run_scan(self,button):
        '''
        Get infotmation and run the simulation
        '''
        pass