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
import gc
import os

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')


class SaveTrajectoryBox:
    """ Class doc """
    
    def __init__ (self, parent):
        """ Class initialiser """
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/trajectory_box.glade'))
        self.builder.connect_signals(self)
        
        self.box = self.builder.get_object('trajectory_box')
        
        self.folder_chooser_button = FolderChooserButton(main =  parent)
        self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
        self.folder_chooser_button.btn.set_sensitive(False)
        
        
        
        
        '''--------------------------------------------------------------------------------------------'''
        self.format_store = Gtk.ListStore(str)
        formats = [
            "pDynamo / pkl" ,
            "amber / crd"               ,
            "charmm / dcd"             ,
            "xyz"   ,
            ]
        for format in formats:
            self.format_store.append([format])
            #print (format)
        self.formats_combo = self.builder.get_object('combobox_format')
        self.formats_combo.set_model(self.format_store)
        #self.formats_combo.connect("changed", self.on_name_combo_changed)
        self.formats_combo.set_model(self.format_store)
        
        renderer_text = Gtk.CellRendererText()
        self.formats_combo.pack_start(renderer_text, True)
        self.formats_combo.add_attribute(renderer_text, "text", 0)
        '''--------------------------------------------------------------------------------------------'''
        self.formats_combo.set_active(0)
        
        
    def on_toggle_save_checkbox (self, widget):
        """ Function doc """
        if self.builder.get_object('checkbox_save_traj').get_active():
            self.builder.get_object('entry_trajectory_name').set_sensitive(True)
            
            self.builder.get_object('label_working_folder').set_sensitive(True)
            #self.builder.get_object('file_chooser_working_folder').set_sensitive(True)
            self.builder.get_object('label_format').set_sensitive(True)
            self.builder.get_object('combobox_format').set_sensitive(True)
            self.builder.get_object('label_trajectory_frequence').set_sensitive(True)
            self.builder.get_object('entry_trajectory_frequence').set_sensitive(True)
            self.folder_chooser_button.btn.set_sensitive(True)

        else:
            self.builder.get_object('entry_trajectory_name').set_sensitive(False)

            self.builder.get_object('label_working_folder').set_sensitive(False)
            #self.builder.get_object('file_chooser_working_folder').set_sensitive(False)
            self.builder.get_object('label_format').set_sensitive(False)
            self.builder.get_object('combobox_format').set_sensitive(False)
            self.builder.get_object('label_trajectory_frequence').set_sensitive(False)
            self.builder.get_object('entry_trajectory_frequence').set_sensitive(False)
            self.folder_chooser_button.btn.set_sensitive(False)
        
    
    def get_folder (self):
        """ Function doc """
        return self.folder_chooser_button.get_folder()
        
    def set_folder (self, folder):
        """ Function doc """
        self.folder_chooser_button.set_folder(folder)


class FolderChooserButton:
    """ Class doc """
    
    def __init__ (self, main = None):
        """ Class initialiser """
        self.main =  main
        self.btn =  Gtk.Button()

        grid = Gtk.Grid ()
        grid.set_column_spacing (10)
        img = Gtk.Image()
        img.set_from_file('easyhybrid/icons/icon_open.png')
        self.label = Gtk.Label ('...')
        self.btn.connect('clicked', self.open_filechooser)



        grid.attach (img, 0, 0, 1, 1)
        grid.attach (self.label, 1, 0, 1, 1)
        grid.show_all ()    

        self.btn.add (grid)        
        #return self.btn 
        self.folder = None
    
    def set_folder (self, folder = '/home'):
        """ Function doc """
        #print('set_folder', folder)
        self.folder = folder
        name = os.path.basename(folder )
        #print( name)
        self.label.set_text(name)
    
    def get_folder (self):
        """ Function doc """
        return self.folder
        
        
    def open_filechooser (self, parent = None):
        """ Function doc """
        dialog = Gtk.FileChooserDialog(
            
            #title="Please choose a file", parent=window, action=Gtk.FileChooserAction.OPEN
            title="Please choose a folder", parent= self.main, action=Gtk.FileChooserAction.SELECT_FOLDER
        )
        dialog.add_buttons(
            Gtk.STOCK_CANCEL,
            Gtk.ResponseType.CANCEL,
            Gtk.STOCK_OPEN,
            Gtk.ResponseType.OK,
        )
        
        if self.folder:
            dialog.set_current_folder(self.folder)
        else:
            dialog.set_current_folder(HOME)
            #self.add_filters(dialog)

        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            #print("Open clicked")
        
            print("File selected: " + dialog.get_filename())
        
            print(dialog.get_filename())
            self.set_folder(folder = dialog.get_filename())
            #print(os.path.dirname( dialog.get_filename() ))

        
        elif response == Gtk.ResponseType.CANCEL:
            print("Cancel clicked")

        dialog.destroy()
        

class GeometryOptimizatrionWindow(Gtk.Window):
    """ Class doc """
    
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/geometry_optimization_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('geometry_optimization_window')
            self.window.set_title('Geometry Optmization Window')
            self.window.set_keep_above(True)
            
            
            '''--------------------------------------------------------------------------------------------'''
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
            self.methods_combo = self.builder.get_object('combobox_geo_opt')
            self.methods_combo.set_model(self.method_store)
            self.methods_combo.connect("changed", self.on_name_combo_changed)
            self.methods_combo.set_model(self.method_store)
            
            renderer_text = Gtk.CellRendererText()
            self.methods_combo.pack_start(renderer_text, True)
            self.methods_combo.add_attribute(renderer_text, "text", 0)
            self.methods_combo.set_active(0)
            '''--------------------------------------------------------------------------------------------'''

            
            
            
 
            
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
            '''--------------------------------------------------------------------------------------------'''

            

            self.save_trajectory_box = SaveTrajectoryBox(parent = self.window)
            self.builder.get_object('geo_opt_parm_box').pack_end(self.save_trajectory_box.box, True, True, 0)
            

            
            # updating data 
            self.update_working_folder_chooser()

            
            
            self.window.show_all()
            self.Visible  = True
    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
    
    #def on_toggle_save_checkbox (self, widget):
    #    """ Function doc """
    #    if self.builder.get_object('checkbox_save_traj').get_active():
    #        self.builder.get_object('entry_trajectory_name').set_sensitive(True)
    #        
    #        self.builder.get_object('label_working_folder').set_sensitive(True)
    #        #self.builder.get_object('file_chooser_working_folder').set_sensitive(True)
    #        self.builder.get_object('label_format').set_sensitive(True)
    #        self.builder.get_object('combobox_format').set_sensitive(True)
    #        self.builder.get_object('label_trajectory_frequence').set_sensitive(True)
    #        self.builder.get_object('entry_trajectory_frequence').set_sensitive(True)
    #        self.folder_chooser_button.btn.set_sensitive(True)
    #
    #    else:
    #        self.builder.get_object('entry_trajectory_name').set_sensitive(False)
    #
    #        self.builder.get_object('label_working_folder').set_sensitive(False)
    #        #self.builder.get_object('file_chooser_working_folder').set_sensitive(False)
    #        self.builder.get_object('label_format').set_sensitive(False)
    #        self.builder.get_object('combobox_format').set_sensitive(False)
    #        self.builder.get_object('label_trajectory_frequence').set_sensitive(False)
    #        self.builder.get_object('entry_trajectory_frequence').set_sensitive(False)
    #        self.folder_chooser_button.btn.set_sensitive(False)

    
    def __init__(self, main = None):
        """ Class initialiser """
        self.easyhybrid_main     = main
        self.Visible             =  False        
        self.residue_liststore = Gtk.ListStore(str, str, str)

        self.opt_methods = { 
                            0 : 'ConjugatedGradient',
                            1 : 'SteepestDescent'   ,
                            2 : 'LFBGS'             ,
                            3 : 'QuasiNewton'       ,
                            4 : 'FIRE'              ,
                             }

    def run_opt (self, button):
        """ Function doc """
    
        '''this combobox has the reference to the starting coordinates of a simulation'''
        combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
        tree_iter = combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
            vismol_object = self.easyhybrid_main.vm_session.vismol_objects_dic[vobject_id]
            
            '''This function imports the coordinates of a vismol_object into the dynamo system in memory.''' 
            self.easyhybrid_main.pDynamo_session.get_coordinates_from_vismol_object_to_pDynamo_system(vismol_object)
        
        
        method_id     = self.builder.get_object('combobox_geo_opt').get_active()
        logFrequency  = int  ( self.builder.get_object('entry_log_frequence').get_text())
        max_int       = int  ( self.builder.get_object('entry_max_int').get_text()  )
        rmsd_tol      = float( self.builder.get_object('entry_rmsd_tol').get_text() )
        
        
        if self.save_trajectory_box.builder.get_object('checkbox_save_traj').get_active():
            folder          = self.save_trajectory_box.folder_chooser_button.get_folder ()
            trajectory_name = self.save_trajectory_box.builder.get_object('entry_trajectory_name').get_text()
            traj_format     = self.save_trajectory_box.builder.get_object('combobox_format').get_active()
            traj_frequence  = int(self.save_trajectory_box.builder.get_object('entry_trajectory_frequence').get_text())
            
            parameters = {
                            
                            'trajectory_name' :  os.path.join(folder, trajectory_name+'.ptGeo')   ,
                            'folder'          :  folder                     ,
                            'optimizer'       :  self.opt_methods[method_id], 
                            'maxIterations'   :  max_int                    ,
                            'log_frequency'   :  logFrequency               ,
                            #'save_pdb'        :  
                            'save_traj'       : True                        , 
                            'save_format'     : '.dcd'                 ,
                            'rmsGradient'     : rmsd_tol                    ,
                            'save_frequency'  : traj_frequence              ,
                            }
            
            '''Simulation routines within easyhybrid receive a parameter dictionary            
            trajectory_name: name to save the trajectory
            maxIterations: maximum number of itetarions (integer) 
            log_frequency: log frequency  (integer)
            save_pdb     : whether to save the final coordinates in pdb format (boolean)
            save_format  : name of the extra binary file ( could be of the format: .dcd, .mdcrd ...) 
            save_frequency : save frame frequency  (integer)
            rmsGradient  : root mean square gradient tolerance ( float )						
            '''            
        else:
            parameters = {
                            
                            'trajectory_name' :  None   ,
                            'folder'          :  None                      ,
                            'optimizer'       :  self.opt_methods[method_id], 
                            'maxIterations'   :  max_int                    ,
                            'log_frequency'   :  logFrequency               ,
                            #'save_pdb'        :  
                            'save_traj'       : False                        , 
                            'save_format'     : None                      ,
                            'rmsGradient'     : rmsd_tol                    ,
                            'save_frequency'  : 0              ,
                            }
            
            
            
            #self.easyhybrid_main.pDynamo_session.run_ConjugateGradientMinimize_SystemGeometry (  
            #                                                                   logFrequency           = logFrequency, 
            #                                                                   maximumIterations      = max_int     ,
            #                                                                   rmsGradientTolerance   = rmsd_tol    , 
            #                                                                   save_trajectory = False,
            #                                                                   trajectory_path = None)
            #
            
            
            
            

        self.easyhybrid_main.pDynamo_session.run_simulation( _parametersList = parameters, 
                                                            _parameters4Plot = None, 
                                                             _simulationType = 'Geometry_Optimization',
                                                             folder          = parameters['folder'])

        self.window.destroy()
        self.Visible    =  False



    
    def on_name_combo_changed(self, widget):
        """ Function doc """
        print('eba - apagar')

    
    def update_working_folder_chooser (self, folder = None):
        """ Function doc """
        if folder:
            #print('update_working_folder_chooser')
            self.save_trajectory_box.set_folder(folder = folder)
        else:
            self.save_trajectory_box.set_folder(folder = HOME)
   
   
   
