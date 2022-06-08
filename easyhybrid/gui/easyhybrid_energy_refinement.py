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



gi.require_version('Gdk', '3.0')
from gi.repository import Gtk, Gdk
from matplotlib.backends.backend_gtk3agg import FigureCanvas  # or gtk3cairo.
from matplotlib.figure import Figure
import numpy as np


from easyhybrid.gui.geometry_optimization_window import SaveTrajectoryBox
from easyhybrid.gui.geometry_optimization_window import FolderChooserButton

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import copy as cp
import gc
import os
from vModel import VismolObject




VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')


class EnergyRefinementWindow():
    
    def OpenWindow (self, vobject = None):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/easyhybrid_energy_refinement.glade'))
            self.builder.connect_signals(self)
            #self.vobject = vobject
            self.window = self.builder.get_object('window')
            self.window.set_title('Energy Refinement Window')
            self.window.connect('destroy', self.CloseWindow)
            
            self.builder.get_object('button_cancel').connect('clicked', self.CloseWindow)
            self.builder.get_object('button_run').connect('clicked', self.on_button_run_clicked)
            
            '''--------------------------------------------------------------------------------------------''' 
            self.combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
            self.starting_coords_liststore = self.main.vm_session.starting_coords_liststore
            self.combobox_starting_coordinates.set_model(self.starting_coords_liststore)
            renderer_text = Gtk.CellRendererText()
            self.combobox_starting_coordinates.pack_start(renderer_text, True)
            self.combobox_starting_coordinates.add_attribute(renderer_text, "text", 0)
            size = len(self.starting_coords_liststore)
            self.combobox_starting_coordinates.set_active(size-1)
            '''--------------------------------------------------------------------------------------------''' 
            
            #'''--------------------------------------------------------------------------------------------------
            self.folder_chooser_button2 = FolderChooserButton(main =  self.window, sel_type = 'file')
            self.builder.get_object('folder_chooser_box2').pack_start(self.folder_chooser_button2.btn, True, True, 0)
            #'''--------------------------------------------------------------------------------------------------
            
            #'''--------------------------------------------------------------------------------------------------
            self.folder_chooser_button1 = FolderChooserButton(main =  self.window)
            self.builder.get_object('folder_chooser_box1').pack_start(self.folder_chooser_button1.btn, True, True, 0)
            #'''--------------------------------------------------------------------------------------------------
            
            
            #'''--------------------------------------------------------------------------------------------------
            self.builder.get_object('box_reaction_coordinate').set_sensitive(False)
            self.builder.get_object('label_CPUs').set_sensitive(False)
            self.builder.get_object('n_CPUs_spinbutton').set_sensitive(False)
            self.builder.get_object('label_input_logfile').set_sensitive(False)
            self.builder.get_object('file_chooser_btn_logfile').set_sensitive(False)
            self.builder.get_object('folder_chooser_box2').set_sensitive(False)
            self.builder.get_object('combobox_coordinate_type').set_sensitive(False)
            self.builder.get_object('label_coordinate_type').set_sensitive(False)
            self.builder.get_object('label_file_or_folder').set_sensitive(False)
            self.builder.get_object('label_input_logfile').set_sensitive(False)
            self.builder.get_object('frame_output').set_sensitive(False)
            #'''--------------------------------------------------------------------------------------------------
            
            #self.window.set_keep_above(True)            
            #self.window.set_default_size(700, 450)
            #
            #self.hbox           = self.builder.get_object('hbox_matplot_figures')
            #self.scale_traj     = self.builder.get_object('scale_trajectory_from_PES')
            #self.adjustment     = Gtk.Adjustment(value         = 0,
            #                                     lower         = 0,
            #                                     upper         = 100,
            #                                     step_increment= 1,
            #                                     page_increment= 1,
            #                                     page_size     = 1)
            #
            #self.scale_traj.set_adjustment ( self.adjustment)
            #self.scale_traj.set_digits(0)
            #self.reaction_coord_label = self.builder.get_object('RC1_RC2_label')
            #
            #
            #
            #
            #
            #'''-------------------------------------------------------------'''
            #self.fig = Figure(figsize=(6, 4))#,constrained_layout=True)
            #self.canvas = FigureCanvas(self.fig)  # a Gtk.DrawingArea
            #self.canvas.mpl_connect('button_press_event', self.on_pick)
            #self.canvas.mpl_connect('motion_notify_event', self.on_motion_notify_event)
            #
            #self.hbox.pack_start(self.canvas, True, True, 0)
            #self.ax = self.fig.add_subplot(1,1,1)
            #'''-------------------------------------------------------------'''
            #
            #
            #
            #
            #'''-------------------------------------------------------------'''
            #self.fig2 = Figure(figsize=(6, 4), dpi=100)
            #self.ax2 = self.fig2.add_subplot(1,1,1)
            #self.ax3 = self.fig2.add_subplot(1,1,1)
            #
            ##self.line2, = self.ax2.plot([], [], '-o')
            ##self.secax = self.ax2.secondary_xaxis('top', functions=(self.seconday_xaxis_function))
            #
            #self.canvas2 = FigureCanvas(self.fig2)  # a Gtk.DrawingArea
            #self.hbox.pack_start(self.canvas2, True, True, 0)        
            #'''-------------------------------------------------------------'''
            #
            #
            #
            #
            #
            #
            #self.grid = self.builder.get_object('grid_setup')
            #
            #sys_selected = 0
            #
            #
            #
            ##'''------------------------------------------------------------------------------------
            #self.vobject_liststore = Gtk.ListStore(str, int)
            #names = [ ]
            #for key , system in self.main.p_session.systems.items():
            #    
            #    for vobject_id in  system['logfile_data'].keys():
            #        print(['vobject_id:', vobject_id ,system['logfile_data'].keys()])
            #        try:
            #            _vobject = self.main.vm_session.vobjects_dic[vobject_id]
            #            print(['_vobject:', _vobject.name,_vobject.index])
            #            self.vobject_liststore.append([_vobject.name, _vobject.index])
            #        except:
            #            print('self.vobject_liststore.append([_vobject.name, _vobject.index])', 'failed!')
            ##self.vobject_liststore.append(['all', _vobject.index])
            #
            #self.coordinates_combobox = Gtk.ComboBox()
            #renderer_text = Gtk.CellRendererText()
            #self.coordinates_combobox.pack_start(renderer_text, True)
            #self.coordinates_combobox.add_attribute(renderer_text, "text", 0)
            #self.coordinates_combobox.set_model(self.vobject_liststore)
            #self.coordinates_combobox.connect('changed', self.on_coordinates_combobox_change)
            #self.grid.attach (self.coordinates_combobox, 1, 0, 1, 1)
            ##------------------------------------------------------------------------------------
            #
            #
            #
            #
            #
            ##------------------------------------------------------------------------------------
            #self.data_combobox = Gtk.ComboBox()
            #renderer_text = Gtk.CellRendererText()
            #self.data_combobox.pack_start(renderer_text, True)
            #self.data_combobox.add_attribute(renderer_text, "text", 0)
            #self.data_combobox.set_model(self.data_liststore)
            #self.data_combobox.connect('changed', self.on_data_combobox_change)
            #self.grid.attach (self.data_combobox, 3, 0, 1, 1)
            ##------------------------------------------------------------------------------------
            #
            #
            #
            ##------------------------------------------------------------------------------------
            #if sys_selected:
            #    self.coordinates_combobox.set_active(sys_selected)
            #else:
            #    self.coordinates_combobox.set_active(0)
            ##------------------------------------------------------------------------------------'''
 
            self.window.show_all()
            self.Visible  = True
    
        else:
            pass
            #_id = self.coordinates_combobox.get_active()
            #self.vobject_liststore = Gtk.ListStore(str, int)
            #names = [ ]
            #for key , system in self.main.p_session.systems.items():
            #    
            #    for vobject_id in  system['logfile_data'].keys():
            #        try:
            #            _vobject = self.main.vm_session.vobjects_dic[vobject_id]
            #            print(['_vobject:', _vobject.name,_vobject.index])
            #            self.vobject_liststore.append([_vobject.name, _vobject.index])
            #        except:
            #            print('self.vobject_liststore.append([_vobject.name, _vobject.index])', 'failed!')
            #
            ##self.vobject_liststore.append(['all', _vobject.index])
            #self.coordinates_combobox.set_model(self.vobject_liststore)
            #self.coordinates_combobox.set_active(_id)


    def on_coordinates_combobox_change (self, widget):
        """ Function doc """
        _id = self.coordinates_combobox.get_active()
        print(_id)
        vobject_index = None
        #-----------------------------------------------------------------------------
        _iter = self.coordinates_combobox.get_active_iter()
        if _iter is not None:
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.coordinates_combobox.get_model()
            _name, vobject_index = model[_iter][:2]
            print ('\n\n\_name, vobject_index:', _name, vobject_index, '\n\n')
        #-----------------------------------------------------------------------------
        self.vobject = self.main.vm_session.vobjects_dic[vobject_index]

        self.data_liststore.clear()
        for index , data in enumerate(self.main.p_session.systems[self.vobject.easyhybrid_system_id]['logfile_data'][vobject_index]):
            #print(data)
            self.data_liststore.append([data['name'], index])
        
        
        #self.data_liststore.append(['all', 2])
        
        self.data_combobox.set_active(0)


    def on_data_combobox_change (self, widget):
        """ Function doc """

        _iter = self.data_combobox.get_active_iter()
        if _iter is not None:
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.data_combobox.get_model()
            _name, index = model[_iter][:2]
            #print ('\n\n\_name, index:', _name,  index, '\n\n')
        
        #self.vobject = self.main.vm_session.vobjects_dic[vobject_index]
        self.data = self.main.p_session.systems[self.vobject.easyhybrid_system_id]['logfile_data'][self.vobject.index][index] 
        print(self.data)
        self._draw_data(cla = True)

    def __init__(self, main = None ):
        """ Class initialiser """
        self.main     = main
        #self.p_session           = self.easyhybrid_main.p_session
        #self.vm_session          = main.vm_session
        self.Visible             =  False        
        
        self.vobject_liststore   = Gtk.ListStore(str, int)
        self.data_liststore      = Gtk.ListStore(str, int)

    
    def radiobutton_single_point_toggled_cb (self, widget):
        """ Function doc """
        if self.builder.get_object('radiobutton_single_point').get_active():
            self.builder.get_object('box_reaction_coordinate').set_sensitive(False)
            self.builder.get_object('label_CPUs').set_sensitive(False)
            self.builder.get_object('n_CPUs_spinbutton').set_sensitive(False)
            self.builder.get_object('label_input_logfile').set_sensitive(False)
            self.builder.get_object('file_chooser_btn_logfile').set_sensitive(False)
            self.builder.get_object('folder_chooser_box2').set_sensitive(False)
            self.builder.get_object('combobox_coordinate_type').set_sensitive(False)
            self.builder.get_object('label_coordinate_type').set_sensitive(False)
            self.builder.get_object('label_file_or_folder').set_sensitive(False)
            self.builder.get_object('label_input_logfile').set_sensitive(False)
            self.builder.get_object('frame_output').set_sensitive(False)
        else:
            self.builder.get_object('box_reaction_coordinate').set_sensitive(True)
            self.builder.get_object('label_CPUs').set_sensitive(True)
            self.builder.get_object('n_CPUs_spinbutton').set_sensitive(True)
            self.builder.get_object('label_input_logfile').set_sensitive(True)
            self.builder.get_object('file_chooser_btn_logfile').set_sensitive(True)
            self.builder.get_object('folder_chooser_box2').set_sensitive(True)
            self.builder.get_object('combobox_coordinate_type').set_sensitive(True)
            self.builder.get_object('label_coordinate_type').set_sensitive(True)
            self.builder.get_object('label_file_or_folder').set_sensitive(True)
            self.builder.get_object('label_input_logfile').set_sensitive(True)
            self.builder.get_object('frame_output').set_sensitive(True)
        
        
    def change_check_button_reaction_coordinate (self, widget):
        """ Function doc """
        #radiobutton_bidimensional = self.builder.get_object('radiobutton_bidimensional')
        if self.builder.get_object('label_check_button_reaction_coordinate2').get_active():
            self.box_reaction_coordinate2.set_sensitive(True)
            self.builder.get_object('n_CPUs_spinbutton').set_sensitive(True)
            self.builder.get_object('n_CPUs_label').set_sensitive(True)
        else:
            self.box_reaction_coordinate2.set_sensitive(False)
            self.builder.get_object('n_CPUs_spinbutton').set_sensitive(False)
            self.builder.get_object('n_CPUs_label')     .set_sensitive(False)

    #======================================================================================
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False


    def on_button_run_clicked (self, button):
        """ Function doc """
        
        #----------------------------------------------------------------------------------------------
        #                            S I N G L E    P O I N T 
        if self.builder.get_object('radiobutton_single_point').get_active():
            
            combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
            tree_iter = combobox_starting_coordinates.get_active_iter()
            if tree_iter is not None:
                
                '''selecting the vismol object from the content that is in the combobox '''
                model = combobox_starting_coordinates.get_model()
                name, vobject_id = model[tree_iter][:2]
                vobject = self.main.vm_session.vobjects_dic[vobject_id]
                
                '''This function imports the coordinates of a vobject into the dynamo system in memory.''' 
                print('vobject:', vobject.name, len(vobject.frames) )
                self.main.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject)

            energy = self.main.p_session.get_energy()
            
            self.CloseWindow ( button =None, data  = None)
            
            dialog = EasyHybridDialogEnergy(parent = self.window, energy = energy, name = vobject.name)
            response = dialog.run()
            dialog.destroy()
        #----------------------------------------------------------------------------------------------








class EasyHybridDialogEnergy(Gtk.Dialog):
    def __init__(self, parent, energy = None, name = 'UNK'):
        super().__init__(title="Energy Dialog", transient_for=parent, flags=0)
        self.add_buttons(
             Gtk.STOCK_OK, Gtk.ResponseType.OK
        )

        self.set_default_size(300, 100)
        string  = "Coordinates = "+name 
        string += "\n\nEnergy = {0:.5f} (KJ/mol)".format(energy)
        
        label = Gtk.Label(label=string)

        box = self.get_content_area()
        box.add(label)
        self.show_all()

