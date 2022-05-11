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
            
            renderer_text = Gtk.CellRendererText()
            self.combobox_starting_coordinates.pack_start(renderer_text, True)
            self.combobox_starting_coordinates.add_attribute(renderer_text, "text", 0)
            
            size = len(self.starting_coords_liststore)
            self.combobox_starting_coordinates.set_active(size-1)

            #'''--------------------------------------------------------------------------------------------'''
            self.md_methods_liststore = Gtk.ListStore(str)
            
            #methods = ["Velocity Verlet Dynamics", "Leap Frog Dynamics","Langevin Dynamics"]
            
            for key, method in self.md_intergators_dict.items():
                self.md_methods_liststore.append([method])
                print (method)
            
            self.methods_combo = self.builder.get_object('comobobox_md_integrator')
            self.methods_combo.set_model(self.md_methods_liststore)
            #self.methods_combo.connect("changed", self.on_name_combo_changed)
            self.methods_combo.set_model(self.md_methods_liststore)
            #
            renderer_text = Gtk.CellRendererText()
            self.methods_combo.pack_start(renderer_text, True)
            self.methods_combo.add_attribute(renderer_text, "text", 0)
            #'''--------------------------------------------------------------------------------------------'''
            self.methods_combo.set_active(0)

            #'''--------------------------------------------------------------------------------------------'''
            self.temp_scale_option_store = Gtk.ListStore(str)
            
            #temp_scale_options = ["constant", "linear","exponential"]
            #temp_scale_options
            for key,temp_scale_option in self.temp_scale_options.items():
                self.temp_scale_option_store.append([temp_scale_option])
                print (temp_scale_option)
            
            self.temp_scale_options_combo = self.builder.get_object('combobox_temperature_scale_options')
            self.temp_scale_options_combo.set_model(self.temp_scale_option_store)
            #self.temp_scale_options_combo.connect("changed", self.on_name_combo_changed)
            self.temp_scale_options_combo.set_model(self.temp_scale_option_store)
            #
            renderer_text = Gtk.CellRendererText()
            self.temp_scale_options_combo.pack_start(renderer_text, True)
            self.temp_scale_options_combo.add_attribute(renderer_text, "text", 0)
            #'''--------------------------------------------------------------------------------------------'''
            self.temp_scale_options_combo.set_active(0)

            #md_parm_box = self.builder.get_object('md_parm_box2')
            
            self.save_trajectory_box = SaveTrajectoryBox(parent = self.window)
            self.builder.get_object('md_parm_box2').pack_end(self.save_trajectory_box.box, True, True, 0)
            self.update_working_folder_chooser()

            #------------------------------------------------------------------------------------
            job_list_canvas = self.builder.get_object('job_list_canvas')
            
            
            '''
            software_list = [
                ("heating"        , '2002', "300"),
                ("Equilibration"  , '2004', "300"),
                ("Data-collection", '2004', "300"),
            ]
                        
            # Creating the ListStore model
            self.job_liststore = Gtk.ListStore(str, str, str)
            for software_ref in software_list:
                self.job_liststore.append(list(software_ref))
            self.current_filter_language = None
            '''
            
            # creating the treeview, making it use the filter as a model, and adding the columns
            self.treeview = Gtk.TreeView(model=self.job_liststore)
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
        self.job_liststore = Gtk.ListStore(str, str, str)
        
        
        self.md_intergators_dict = {
                                   0:"Velocity Verlet Dynamics", 
                                   1:"Leap Frog Dynamics",
                                   2:"Langevin Dynamics"
                                   }
                       
        self.temp_scale_options = {0: "constant", 
                                   1: "linear",
                                   2: "exponential"}
        
        
        self.joblist = []
        #self.opt_methods = { 
        #                    0 : 'ConjugatedGradient',
        #                    1 : 'SteepestDescent'   ,
        #                    2 : 'LFBGS'             ,
        #                    3 : 'QuasiNewton'       ,
        #                    4 : 'FIRE'              ,
        #                     }




    def add_job_to_joblist (self, button):
        '''
        Set up and execute molecular dynamics simulations.:
        Mandatory keys in self.parameters:
            "MD_method"	 	 : string containing the integrator algorithm name
            "protocol" 		 : string indicating if is a normal run or for heating
            "nsteps"   		 : Number of steps to be taken in the simulation
            "trajectory_name":
        Optinal  :
            "temperature" 			  : float with the simulation temperature. If not passed we assume 300.15K as default.
            "coll_freq"  			  : integer with the colision frequency. Generally set for Langevin integrator. 
            "pressure"   			  : float with the simulation pressure. If not passed we assume 1.0bar as default.
            "pressure_coupling"		  : boolean indicating if is to control the simulation pressure.
            "temperature_scale_option": string with the type of temperature scaling. Default is 'linear' ( relevant for "heating" protocol)
            "temperature_scale"		  : float with the  temperature scaling step. Default is 10K  ( relevant for "heating" protocol)
            "start_temperatue"		  : float with the start temperature for heating protocol
            "timeStep"   			  : float indicating the size of integration time step. 0.001 ps is taken as default.					
            "sampling_factor"		  : integer indicating in which frequency to save/collect structure/data. default 0.
            "seed"					  : integer indicating the seed for rumdomness of the simulations.
            "log_frequency"     	  : integer indicating the frequency of the screen log output.
        plot parameters keys in self.parameters
            Optinal   :
            "show"					: whether to show the analysis plots in the simulation end.
            "calculate_distances"	: indicate if to calculate distances distributions of passed reaction coordinates
            "ATOMS_RC1"             : list of atoms for the first reaction coordinate to be analyzed 
            "ATOMS_RC2"             : list of atoms for the second reaction coordinate to be analyzed 
        '''		
        
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
            vismol_object = self.easyhybrid_main.vm_session.vismol_objects_dic[vobject_id]
            
            '''This function imports the coordinates of a vismol_object into the dynamo system in memory.''' 
            print('vismol_object:', vismol_object.name, len(vismol_object.frames) )
            self.easyhybrid_main.p_session.get_coordinates_from_vismol_object_to_pDynamo_system(vismol_object)


        coord_id            = self.builder.get_object('combobox_starting_coordinates').get_active()
        integrator_id       = self.builder.get_object('comobobox_md_integrator').get_active()
        number_of_steps     = int(self.builder.get_object('entry_number_of_steps').get_text())
        temp_scale_id       = self.builder.get_object('combobox_temperature_scale_options').get_active()
        temp_start          = int(self.builder.get_object('entry_temp_start').get_text())
        temp_end            = int(self.builder.get_object('entry_temp_end').get_text())
        temp_scale_factor   = int(self.builder.get_object('entry_temp_scale_factor').get_text())
        time_step           = float(self.builder.get_object('entry_time_step').get_text())
        log_frequence       = int(self.builder.get_object('entry_log_frequence').get_text())
        random_seed         = int(self.builder.get_object('entry_random_seed').get_text())
        collision_frequency = float(self.builder.get_object('entry_collision_frequency').get_text())
        
        
        if self.builder.get_object('check_pressure_control').get_active():
            pressure_control = True
        else:
            pressure_control = False
            
        
        pressure          = float(self.builder.get_object('entry_pressure').get_text())
        pressure_coupling = float(self.builder.get_object('entry_pressure_coupling').get_text())
        
        parameters = {
                    'coord_id'           : coord_id,           
                    'integrator_id'      : integrator_id,        
                    'number_of_steps'    : number_of_steps,      
                    'temp_scale_id'      : temp_scale_id,        
                    'temp_start'         : temp_start,           
                    'temp_end'           : temp_end,             
                    'temp_scale_factor'  : temp_scale_factor,    
                    'time_step'          : time_step,            
                    'log_frequence'      : log_frequence,       
                    'random_seed'        : random_seed,          
                    'collision_frequency': collision_frequency,  
                    'pressure_control'   : pressure_control,   
                    'pressure'           : pressure,          
                    'pressure_coupling'  : pressure_coupling, 
                    }
        #print(parameters)
        
        self.joblist.append(parameters)
        
        
        self.job_liststore.append.clear()
        
        for job in job_liststore:
            self.job_liststore.append(list([job, job ,job]))
        #self.job_liststore.append(list(['aqui', '123' ,'cocozao']))




    def update_working_folder_chooser (self, folder = None):
        """ Function doc """
        if folder:
            print('update_working_folder_chooser')
            self.save_trajectory_box.set_folder(folder = folder)
        else:
            self.save_trajectory_box.set_folder(folder = HOME)

    def combobox_change (self, widget):
        """ Function doc """
        if widget  == self.builder.get_object('comobobox_md_integrator'):
            
            #velocity verlet molecular dynamics
            if 0 == widget.get_active():
                self.builder.get_object('check_pressure_control').set_sensitive(False)
                self.builder.get_object('entry_pressure').set_sensitive(False)
                
                self.builder.get_object('label_pressure_coupling').set_sensitive(False)
                self.builder.get_object('entry_pressure_coupling').set_sensitive(False)
                
                self.builder.get_object('label_collision_freq').set_sensitive(False)
                self.builder.get_object('entry_collision_frequency').set_sensitive(False)
            
            #Leap Frog molecular dynamics
            elif 1 == widget.get_active():
                self.builder.get_object('check_pressure_control').set_sensitive(True)
                self.builder.get_object('entry_pressure').set_sensitive(True)
                
                self.builder.get_object('label_pressure_coupling').set_sensitive(True)
                self.builder.get_object('entry_pressure_coupling').set_sensitive(True)
                
                self.builder.get_object('label_collision_freq').set_sensitive(False)
                self.builder.get_object('entry_collision_frequency').set_sensitive(False)
            
            #Langevin molecular dynamics
            elif 2 == widget.get_active():
                self.builder.get_object('check_pressure_control').set_sensitive(False)
                self.builder.get_object('entry_pressure').set_sensitive(False)
                
                self.builder.get_object('label_pressure_coupling').set_sensitive(False)
                self.builder.get_object('entry_pressure_coupling').set_sensitive(False)
                
                self.builder.get_object('label_collision_freq').set_sensitive(True)
                self.builder.get_object('entry_collision_frequency').set_sensitive(True)
            else:
                pass

        
        
        
        elif widget == self.builder.get_object('combobox_temperature_scale_options'):
            temp_scale_id = widget.get_active()
            
            if temp_scale_id == 0:
                self.builder.get_object('entry_temp_end').set_sensitive(False)
                self.builder.get_object('label_temp_end').set_sensitive(False)
            else:
                self.builder.get_object('entry_temp_end').set_sensitive(True)
                self.builder.get_object('label_temp_end').set_sensitive(True)

        else:
            pass
