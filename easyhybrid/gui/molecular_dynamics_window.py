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
from easyhybrid.gui.geometry_optimization_window import SaveTrajectoryBox
from easyhybrid.gui.geometry_optimization_window import FolderChooserButton
import gc
import os

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')


class MolecularDynamicsWindow():
    
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/easyhybrid_molecular_dynamics_window.glade'))
            self.builder.connect_signals(self)
            #
            self.window = self.builder.get_object('molecular_dynamics_window')
            self.window.set_title('Molecular Dynamics Window')
            #self.window.set_keep_above(True)

            '''--------------------------------------------------------------------------------------------'''
            self.combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
            self.starting_coords_liststore = self.easyhybrid_main.vm_session.starting_coords_liststore
            self.combobox_starting_coordinates.set_model(self.starting_coords_liststore)
            #self.combobox_starting_coordinates.connect("changed", self.on_name_combo_changed)
            self.combobox_starting_coordinates.set_model(self.starting_coords_liststore)
            
            renderer_text = Gtk.CellRendererText()
            self.combobox_starting_coordinates.pack_start(renderer_text, True)
            self.combobox_starting_coordinates.add_attribute(renderer_text, "text", 0)
            
            size = len(self.starting_coords_liststore)
            self.combobox_starting_coordinates.set_active(size-1)

            #'''--------------------------------------------------------------------------------------------'''
            self.method_store = Gtk.ListStore(str)
            
            methods = ["Velocity Verlet Dynamics", "Leap Frog Dynamics","Langevin Dynamics"]
            
            for method in methods:
                self.method_store.append([method])
                print (method)
            
            self.methods_combo = self.builder.get_object('md_integrator_comobobox')
            self.methods_combo.set_model(self.method_store)
            #self.methods_combo.connect("changed", self.on_name_combo_changed)
            self.methods_combo.set_model(self.method_store)
            #
            renderer_text = Gtk.CellRendererText()
            self.methods_combo.pack_start(renderer_text, True)
            self.methods_combo.add_attribute(renderer_text, "text", 0)
            #'''--------------------------------------------------------------------------------------------'''
            self.methods_combo.set_active(0)

            #'''--------------------------------------------------------------------------------------------'''
            self.temp_scale_option_store = Gtk.ListStore(str)
            
            temp_scale_options = ["constant", "linear","exponential"]
            
            for temp_scale_option in temp_scale_options:
                self.temp_scale_option_store.append([temp_scale_option])
                print (temp_scale_option)
            
            self.temp_scale_options_combo = self.builder.get_object('temperature_scale_option_combobox')
            self.temp_scale_options_combo.set_model(self.temp_scale_option_store)
            #self.temp_scale_options_combo.connect("changed", self.on_name_combo_changed)
            self.temp_scale_options_combo.set_model(self.temp_scale_option_store)
            #
            renderer_text = Gtk.CellRendererText()
            self.temp_scale_options_combo.pack_start(renderer_text, True)
            self.temp_scale_options_combo.add_attribute(renderer_text, "text", 0)
            #'''--------------------------------------------------------------------------------------------'''
            self.temp_scale_options_combo.set_active(0)

            md_parm_box = self.builder.get_object('md_parm_box')
            
            self.save_trajectory_box = SaveTrajectoryBox(parent = self.window)
            self.builder.get_object('md_parm_box').pack_end(self.save_trajectory_box.box, True, True, 0)
            self.update_working_folder_chooser()

            #------------------------------------------------------------------------------------
            job_list_canvas = self.builder.get_object('job_list_canvas')
            
            software_list = [
                ("heating"        , '2002', "300"),
                ("Equilibration"  , '2004', "300"),
                ("Data-collection", '2004', "300"),
            ]
                        
            # Creating the ListStore model
            self.software_liststore = Gtk.ListStore(str, str, str)
            for software_ref in software_list:
                self.software_liststore.append(list(software_ref))
            self.current_filter_language = None
            
            # creating the treeview, making it use the filter as a model, and adding the columns
            self.treeview = Gtk.TreeView(model=self.software_liststore)
            for i, column_title in enumerate(
                ["Job", "nSteps", "temp"]
            ):
                renderer = Gtk.CellRendererText()
                column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                self.treeview.append_column(column)

            self.scrollable_treelist = Gtk.ScrolledWindow()
            self.scrollable_treelist.set_vexpand(True)
            #self.grid.attach(self.scrollable_treelist, 0, 0, 8, 10)
            #self.grid.attach_next_to(
            #    self.buttons[0], self.scrollable_treelist, Gtk.PositionType.BOTTOM, 1, 1
            #)
            #for i, button in enumerate(self.buttons[1:]):
            #    self.grid.attach_next_to(
            #        button, self.buttons[i], Gtk.PositionType.RIGHT, 1, 1
            #    )
            self.scrollable_treelist.add(self.treeview)
            
            job_list_canvas.add(self.scrollable_treelist)

            self.window.show_all()
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
        
        #self.opt_methods = { 
        #                    0 : 'ConjugatedGradient',
        #                    1 : 'SteepestDescent'   ,
        #                    2 : 'LFBGS'             ,
        #                    3 : 'QuasiNewton'       ,
        #                    4 : 'FIRE'              ,
        #                     }

    def add_job_to_list (self, button):
        """ Function doc """
        self.software_liststore.append(list(['aqui', '123' ,'cocozao']))

    def update_working_folder_chooser (self, folder = None):
        """ Function doc """
        if folder:
            print('update_working_folder_chooser')
            self.save_trajectory_box.set_folder(folder = folder)
        else:
            self.save_trajectory_box.set_folder(folder = HOME)
