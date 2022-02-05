#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  __main__.py
#  
#  Copyright 2016 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
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

import gi, sys
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk
#import os

#w = Gtk.Window()
#f = Gtk.Image()
#f.set_from_file("/home/fernando/Pictures/Screenshot from 2021-06-12 09-49-19.png")
#w.add(f)
#w.show_all()
#Gtk.main()

from vCore.VismolSession  import VisMolSession
#from GTKGUI               import VismolMain 
from easyhybrid.GUI       import EasyHybridMainWindow
#from easyhybrid.GUI       import LabelWindow
from vModel.Representations   import LinesRepresentation
from vModel.Representations   import NonBondedRepresentation
from vModel.Representations   import SticksRepresentation
from vModel.Representations   import DotsRepresentation
from vModel.Representations   import SpheresRepresentation
from vModel.Representations   import GlumpyRepresentation
from vModel.Representations   import RibbonsRepresentation
from vModel.Representations   import SurfaceRepresentation
from vModel.Representations   import WiresRepresentation
from vModel.Representations   import LabelRepresentation

from vModel.Representations   import DynamicBonds
from vModel.Representations   import CartoonRepresentation

from easyhybrid.Serialization import LoadAndSaveFiles

from vModel import VismolObject
import pickle



class EasyHybridVismolSession(VisMolSession, LoadAndSaveFiles):
    """ Class doc """
   
    #ef __init__ (self, glwidget = False, toolkit = 'gtk3', main_session = None):
    #   """ Function doc """
    #   super().__init__( glwidget = False, toolkit = 'gtk3', main_session = None)
    #   self.treestore = Gtk.TreeStore(
    #                                   str , # Name
    #                                   
    #                                   bool, # toggle active=1
    #                                   bool, # radio  active=2
    #                                   
    #                                   bool, # toggle visible = 3
    #                                   bool, # radio  visible = 4
    #                                   
    #                                   int , # vismol_object index
    #                                   bool, # is vismol_object index visible?
    #                                   int , # pdynamo system index 
    #                                   bool, # is pdynamo system index visible?
    #
    #                                   
    #                                   )
    #   
    def save_serialization_file (self, filename = 'session.easy'):
        """ Function doc """
        #serialization = LoadAndSaveFiles(self, self.main_session.pDynamo_session)
        self.save_session(filename)


    def build_index_list_from_atom_selection (self):
        """ Function doc """
        selection         = self.selections[self.current_selection]
        
        index_list = []                
        for atom in selection.selected_atoms:
            #print(atom.Vobject.easyhybrid_system_id , pdmsys_active)
            true_or_false = self.check_selected_atom (atom)
            if true_or_false:
                index_list.append(atom.index -1)
            else:
                return False
        return index_list


    def check_selected_atom(self, atom, dialog = True):
        '''checks if selected atoms belong to the dynamo system in memory'''
        if atom.Vobject.easyhybrid_system_id != self.main_session.pDynamo_session.active_id:
            #print(atom.index-1, atom.name, atom.resn)
            
            name = self.main_session.pDynamo_session.systems[self.main_session.pDynamo_session.active_id]['name']
            
            dialog = Gtk.MessageDialog(
                        transient_for = self.main_session.window,
                        flags=0,
                        message_type=Gtk.MessageType.INFO,
                        buttons=Gtk.ButtonsType.OK,
                        text="Invalid Atom Selection",
                        )
            dialog.format_secondary_text(
"""Your atom selection does not belong to the active pDynamo system:
 
{} ({}) 

You can choose the active pDynamo system by changing the radio 
button position in the main treeview (active column).""".format(name,self.main_session.pDynamo_session.active_id)
            )
            dialog.run()
            print("INFO dialog closed")
            dialog.destroy()
            return False
        else:
            return True
        
    
    def load_easyhybrid_serialization_file (self, filename):
        """ Function doc """
        #new_session = self.restart_session(filename)
        #serialization = LoadAndSaveFiles(self, self.main_session.pDynamo_session)
        '''
        #--------------------------------------------------------------------------
        self.vismol_objects     = [] # old Vobjects - include molecules
        self.vismol_objects_dic = {} # old Vobjects dic - include molecules
        self.vobj_counter       = 0  # Each vismol object has a unique access key (int), which is generated in the method: add_vismol_object_to_vismol_session.
        self.vismol_geometric_object     = []
        
        self.vismol_geometric_object_dic = {
                                           'pk1pk2' :  None,
                                           'pk2pk3' :  None,
                                           'pk3pk4' :  None,
                                           }
        
        self.atom_id_counter  = 0  # 
        self.atom_dic_id      = {
                                # atom_id : obj_atom 
                                 }
        
        #--------------------------------------------------------------------------
        '''

        
        #print('loading easyhybrid session')
        self.load_session(filename)

    
    
    
    
    def add_vismol_object_to_vismol_session (self, pdynamo_session = None, 
                                                   rep             = {'lines': [], 'nonbonded': []}, 
                                                   vismol_object   = None, 
                                                   vobj_count      = True,
                                                   autocenter      = True):
        """ Function doc """
       
        if vobj_count:
            vismol_object.index = self.vobj_counter
        else:
            pass
        
        #self.vismol_objects.append(vismol_object)
        self.vismol_objects_dic[vismol_object.index] = vismol_object
        #self.append_vismol_object_to_vismol_objects_listStore(vismol_object)
        
        if vobj_count:
            self.vobj_counter += 1
        
            
        vobj_index = vismol_object.index
        sys_index  = vismol_object.easyhybrid_system_id
        
                
        
        if sys_index in self.parents.keys():
            pass
        else:          
            # Creates a new "parent" when a new system is loaded into memory. 
            for row in self.treestore:
                #row[2] = row.path == selected_path
                row[3] =  False
                #row[5] = False 
                for i,j in enumerate(row):
                    print(i, j,)
            
                    
            
            '''
            self.parents[sys_index] = self.treestore.append(None,                                                
                                                           [pdynamo_session.systems[sys_index]['name'], # Name
                                                            False,                                      # toggle active=1
                                                            True,                                       # radio  active=2
                                                            False,                                      # toggle visible = 3
                                                            True,                                       # radio  visible = 4
                                                            False,                                      # radio  active = 5
                                                            False,                                      # is trajectory radio visible?
                                                            vismol_object.easyhybrid_system_id,         # pdynamo system index
                                                            False])                                     # is pdynamo system index visible?
            '''
            
            
            self.parents[sys_index] = self.treestore.append(None,                                                
                                                           
                                                           [pdynamo_session.systems[sys_index]['name'], # Name
                                                            False,                                      # toggle active=1
                                                            False,                                      # toggle visible = 3
                                                            
                                                            True ,                                      # radio  active  = 2
                                                            True ,                                      # radio  visible = 4

                                                            False,                                      # traj radio  active = 5
                                                            False,                                      # is trajectory radio visible?
                                                            
                                                            vismol_object.index,                        #
                                                            vismol_object.easyhybrid_system_id,         # pdynamo system index
                                                            0])                                     # is pdynamo system index visible?
            
            self.gtk_treeview_iters.append(self.parents[sys_index])
            
            
        #or row in self.treestore:
        #   #row[2] = row.path == selected_path
        #   row[5] = False 
        #   for i,j in enumerate(row):
        #       print(i, j,)
            
        for treeview_iter in self.gtk_treeview_iters:
            self.treestore[treeview_iter][5] = False
            #print(self.treestore[treeview_iter][0])
            
        treeview_iter = self.treestore.append(self.parents[vismol_object.easyhybrid_system_id]      ,        #parent
                                          
                                          [vismol_object.name, 
                                           vismol_object.active ,   # toggle active   =1       
                                           True ,                   # toggle visible  = 2                  
                                           
                                           False ,                  # radio  active  = 3                       
                                           False ,                  # radio  visible = 4                      
                                           
                                           True  ,                  # traj radio  active = 5                     
                                           True  ,                  # is trajectory radio visible?  6                   
                                           
                                           vismol_object.index,     # 7
                                           vismol_object.easyhybrid_system_id,   # pdynamo system index  8    
                                           len(vismol_object.frames)] # is pdynamo system index visible?  9 
                                            )
        self.gtk_treeview_iters.append(treeview_iter)
        
        #else:
        #    self.treestore.append(self.parents[vismol_object.easyhybrid_system_id], [vismol_object.name, vismol_object.active , False , False, False, vobj_index, True, vismol_object.easyhybrid_system_id, False])


        print('\n\n\n')

        for vobj in self.vismol_objects:
            print(vobj.name, vobj.easyhybrid_system_id)
        print('\n\n\n')
        print(vismol_object.index)
        
        '''
        if rep:
            vismol_object.create_new_representation (rtype = 'lines')
            self.vismol_objects[-1].create_new_representation(rtype = 'dotted_lines')

            rep  = NonBondedRepresentation (name = 'nonbonded', active = True, _type = 'mol', visObj = vismol_object, glCore = self.glwidget.vm_widget)
            vismol_object.representations[rep.name] = rep
        '''
            

        for key in rep.keys():
            
            if key == 'lines':
                if rep[key] == []:
                    vismol_object.create_new_representation (rtype = 'lines')
                else:
                    
                    vismol_object.create_new_representation (rtype = 'lines', indexes = rep[key])
           
            if key == 'nonbonded':
                if rep[key] == []:
                    vismol_object.create_new_representation (rtype = 'nonbonded')
                else:
                    vismol_object.create_new_representation (rtype = 'nonbonded', indexes = rep[key])
                
            if key == 'sticks':
                if rep[key] == []:
                    vismol_object.create_new_representation (rtype = 'sticks')
                else:
                    vismol_object.create_new_representation (rtype = 'sticks', indexes = rep[key])
            
            if key == 'spheres':
                if rep[key] == []:
                    vismol_object.create_new_representation (rtype = 'spheres')
                else:
                    print(key,rep[key])
                    vismol_object.create_new_representation (rtype = 'spheres', indexes = rep[key])
                
            
                
                
                #self.vismol_objects[-1].generate_indexesresentations (reps_list = self.indexes)
                #print (self.vismol_objects[-1].representations)

            #rep =  CartoonRepresentation(name = 'cartoon', active = True, _type = 'mol', visObj = vismol_object, glCore = self.glwidget.vm_widget)
            #vismol_object.representations[rep.name] = rep
            
            
            #vismol_object.create_new_representation (rtype = 'ribbons')
            

        if autocenter:
            #print(self.vismol_objects[-1].mass_center)
            self.glwidget.vm_widget.center_on_coordinates(vismol_object, vismol_object.mass_center)
        else:
            self.glwidget.vm_widget.queue_draw()
        self.gtk_widgets_update ()

        
        print('vismol_objects_dic', self.vismol_objects_dic)




    
    def insert_glmenu (self, bg_menu  = None, 
                            sele_menu = None, 
                             obj_menu = None, 
                            pick_menu = None):
        """ Function doc """
        



        def _viewing_selection_mode_atom (_):
            """ Function doc """
            self.viewing_selection_mode(sel_type = 'atom')
        def _viewing_selection_mode_residue (_):
            """ Function doc """
            self.viewing_selection_mode(sel_type = 'residue')
        def _viewing_selection_mode_chain (_):
            """ Function doc """
            self.viewing_selection_mode(sel_type = 'chain')

        def _selection_type_picking(_):
            
            if self.selection_box_frane:
                self.selection_box_frane.change_toggle_button_selecting_mode_status(True)
            else:
                self._picking_selection_mode = True
            self.glwidget.queue_draw()
        
        def _selection_type_viewing(_):
            if self.selection_box_frane:
                self.selection_box_frane.change_toggle_button_selecting_mode_status(False)
            else:
                self._picking_selection_mode = False
            self.glwidget.queue_draw()

        if sele_menu is None:
            ''' Standard Sele Menu '''
            
            def dynamic_test (_):
                """ Function doc """
                self.teste3()
            
            def select_test (_):
                """ Function doc """
                self.select(indexes = 'all')
            
            def menu_show_lines (_):
                """ Function doc """
                self.show_or_hide( _type = 'lines', show = True)

            def menu_hide_lines (_):
                """ Function doc """
                #print('hide')
                self.show_or_hide( _type = 'lines', show = False)

            def menu_show_sticks (_):
                """ Function doc """
                self.show_or_hide( _type = 'sticks', show = True)
            
            def menu_show_nonbonded (_):
                """ Function doc """
                self.show_or_hide( _type = 'nonbonded', show = True)
            
            def menu_hide_nonbonded (_):
                """ Function doc """
                self.show_or_hide( _type = 'nonbonded', show = False)

            def menu_hide_sticks (_):
                """ Function doc """
                self.show_or_hide( _type = 'sticks', show = False)

            def menu_show_spheres (_):
                """ Function doc """
                self.show_or_hide( _type = 'spheres', show = True)

            def menu_hide_spheres (_):
                """ Function doc """
                self.show_or_hide( _type = 'spheres', show = False)
            
            def menu_show_dots (_):
                """ Function doc """
                self.show_or_hide( _type = 'dots', show = True)

            def menu_hide_dots (_):
                """ Function doc """
                self.show_or_hide( _type = 'dots', show = False)
            
            def set_as_qc_atoms (_):
                """ Function doc """
                #selection = self.selections[self.current_selection]
                pdmsys_active = self.main_session.pDynamo_session.active_id
                qc_list = self.build_index_list_from_atom_selection()
                if qc_list:
                    self.main_session.pDynamo_session.systems[pdmsys_active]['qc_table'] = qc_list
                    self.main_session.run_dialog_set_QC_atoms()

            def set_as_free_atoms (_):
                """ Function doc """
                selection         = self.selections[self.current_selection]
                
                freelist = []                
                for atom in selection.selected_atoms:
                    #print(atom.index, atom.name, atom.color) 
                    true_or_false = self.check_selected_atom(atom, dialog = True)
                    
                    if true_or_false:
                        freelist.append(atom.index -1)
                        atom.color = atom.init_color(atom.symbol) 
                        #true_or_false = self.check_selected_atom( atom, dialog = True)
                    else:
                        return False
                    
                    #print(atom.index, atom.name, atom.color) 
                #----------------------------------------------
                pdmsys_active =   self.main_session.pDynamo_session.active_id
                #fixedlist = fixedlist + self.main_session.pDynamo_session.systems[pdmsys_active]['fixed_table']
                a = set(self.main_session.pDynamo_session.systems[pdmsys_active]['fixed_table'])
                b = set(freelist)
                
                c = a - b
                
                print (a)
                print (b)
                #Combining with list that the already exists  
                fixedlist =  set(self.main_session.pDynamo_session.systems[pdmsys_active]['fixed_table']) -set(freelist)
                #guarantee that the atom index appears only once in the list
                fixedlist = list(c) 
                print ('fixedlist',fixedlist)
                #sending to pDynamo
                refresh = self.main_session.pDynamo_session.define_free_or_fixed_atoms_from_iterable (fixedlist)
                if refresh:
                    #self.main_session.pDynamo_session.refresh_qc_and_fixed_representations()
                    self.glwidget.vm_widget.queue_draw()
                #self.main_session.pDynamo_session.vismol_selection_qc = selection.copy()
            
            def prune_atoms (_):
                """ Function doc """
                fixedlist = self.build_index_list_from_atom_selection()
                if fixedlist:
                    fixedlist = list(set(fixedlist))
                    #self.main_session.pDynamo_session.define_free_or_fixed_atoms_from_iterable (fixedlist)
                    self.main_session.pDynamo_session.prune_system (selection = fixedlist, label = 'Pruned System', summary = True)
            
            def set_as_fixed_atoms (_):
                """ Function doc """
                
                fixedlist = self.build_index_list_from_atom_selection()
                
                if fixedlist:
                    pdmsys_active = self.main_session.pDynamo_session.active_id
                    fixedlist = fixedlist + self.main_session.pDynamo_session.systems[pdmsys_active]['fixed_table']
                    #guarantee that the atom index appears only once in the list
                    fixedlist = list(set(fixedlist)) 
                    print ('fixedlist',fixedlist)
                    #sending to pDynamo
                    refresh = self.main_session.pDynamo_session.define_free_or_fixed_atoms_from_iterable (fixedlist)
                    if refresh:
                        self.glwidget.vm_widget.queue_draw()
                    #self.main_session.pDynamo_session.vismol_selection_qc = selection.copy()
            
            def invert_selection (_):
                """ Function doc """
                #print('self.selections[self.current_selection].invert_selection()')
                self.selections[self.current_selection].invert_selection()
            
            
            sele_menu = { 
                    'header' : ['MenuItem', None],
                    
                    
                    
                    'separator1':['separator', None],
                    
                    
                    'show'   : [
                                'submenu' ,{
                                            
                                            'lines'         : ['MenuItem', menu_show_lines],
                                            'sticks'        : ['MenuItem', menu_show_sticks],
                                            'spheres'       : ['MenuItem', menu_show_spheres],
                                            'dots'          : ['MenuItem', menu_show_dots],
                                            'dynamic bonds' : ['MenuItem', dynamic_test],
                                            'separator2'    : ['separator', None],
                                            'nonbonded'     : ['MenuItem', menu_show_nonbonded],
                    
                                           }
                               ],
                    
                    
                    'hide'   : [
                                'submenu',  {
                                            'lines'    : ['MenuItem', menu_hide_lines],
                                            'sticks'   : ['MenuItem', menu_hide_sticks],
                                            'spheres'  : ['MenuItem', menu_hide_spheres],
                                            'dots'     : ['MenuItem', menu_hide_dots],
                                            'separator2'    : ['separator', None],
                                            'nonbonded': ['MenuItem', menu_hide_nonbonded],
                                            }
                                ],
                    
                    'Invert Selection':['MenuItem', invert_selection],
                    
                    'separator2':['separator', None],

                    
                    
                    'Selection type'   : [
                                'submenu' ,{
                                            
                                            'viewing'   :  ['MenuItem', _selection_type_viewing],
                                            'picking'   :  ['MenuItem', _selection_type_picking],
                                            #'separator2':['separator', None],
                                            #'nonbonded' : ['MenuItem', None],
                    
                                           }
                                        ],
                    
                    'Selection Mode'   : [
                                'submenu' ,{
                                            
                                            'Atoms'     :  ['MenuItem', _viewing_selection_mode_atom],
                                            'Residue'   :  ['MenuItem', _viewing_selection_mode_residue],
                                            'Chain'     :  ['MenuItem', _viewing_selection_mode_chain],
                                            #'separator2':['separator', None],
                                            #'nonbonded' : ['MenuItem', None],
                    
                                           }
                               ],
                    
                    'separator3':['separator', None],
                    
                    'Set as QC atoms2'      :  ['MenuItem', set_as_qc_atoms],
                    
                    'separator4':['separator', None],

                    'Set as fixed atoms'   :  ['MenuItem', set_as_fixed_atoms],
                    'Set as free atoms'   :  ['MenuItem', set_as_free_atoms],
                    
                    'separator5':['separator', None],
                    'prune to selection'  :  ['MenuItem', prune_atoms],

                    'separator6':['separator', None],

                    
                    'Label Mode':  ['submenu' , {
                                            'Atom'         : [
                                                               'submenu', {
                                                                           'lines'    : ['MenuItem', None],
                                                                           'sticks'   : ['MenuItem', None],
                                                                           'spheres'  : ['MenuItem', None],
                                                                           'nonbonded': ['MenuItem', None],
                                                                           }
                                                              ],
                                            
                                            'Atom index'   : ['MenuItem', None],
                                            'residue name' : ['MenuItem', None],
                                            'residue_index': ['MenuItem', None],
                                           },
                               ]
                    }
      
        if bg_menu is None:
            ''' Standard Bg Menu'''
            
            def open_structure_data (_):
                """ Function doc """
                #print('ebaaaa')
                self.filechooser   = FileChooser()
                filename = self.filechooser.open()
                self.load (filename, widget = None, autocenter = True)


                
            bg_menu = { 
                    'separator0'   :['separator', None],

                    'Open File'    : ['MenuItem', open_structure_data],
                    
                    'select' : ['MenuItem', select_test],

                    'funcao teste' : ['MenuItem', self.teste],                  
                    'funcao teste2': ['MenuItem', self.teste2], 

                    'separator1':['separator', None],


                    'Selection type'   : [
                                'submenu' ,{
                                            
                                            'viewing'   :  ['MenuItem', _selection_type_viewing],
                                            'picking'   :  ['MenuItem', _selection_type_picking],
                                            #'separator2':['separator', None],
                                            #'nonbonded' : ['MenuItem', None],
                    
                                           }
                                        ],
                    
                    'Selection Mode'   : [
                                'submenu' ,{
                                            
                                            'atoms'     :  ['MenuItem', _viewing_selection_mode_atom],
                                            'residue'   :  ['MenuItem', _viewing_selection_mode_residue],
                                            'chain'     :  ['MenuItem', _viewing_selection_mode_chain],
                                            #'separator2':['separator', None],
                                            #'nonbonded' : ['MenuItem', None],
                    
                                           }
                               ],
                    
                    
                    'hide'   : [
                                'submenu',  {
                                            'lines'    : ['MenuItem', menu_hide_lines],
                                            'sticks'   : ['MenuItem', menu_hide_sticks],
                                            'spheres'  : ['MenuItem', menu_hide_spheres],
                                            'nonbonded': ['MenuItem', None],
                                            }
                                ],
                    
                    
                    'separator2':['separator', None],

                    
                    
                    'label':  ['submenu' , {
                                            'Atom'         : [
                                                               'submenu', {
                                                                           'lines'    : ['MenuItem', None],
                                                                           'sticks'   : ['MenuItem', None],
                                                                           'spheres'  : ['MenuItem', None],
                                                                           'nonbonded': ['MenuItem', None],
                                                                           }
                                                              ],
                                            
                                            'Atom index'   : ['MenuItem', None],
                                            'residue name' : ['MenuItem', None],
                                            'residue_index': ['MenuItem', None],
                                           },
                               ]
                    }

        if obj_menu is None:
            ''' Standard Obj Menu'''
            obj_menu = { 
                    'OBJ menu' : ['MenuItem', None],
                    
                    
                    'separator1':['separator', None],
                    
                    
                    'show'   : [
                                'submenu' ,{
                                            
                                            'lines'    : ['MenuItem', menu_show_lines],
                                            'sticks'   : ['MenuItem', menu_show_sticks],
                                            'spheres'  : ['MenuItem', menu_show_spheres],
                                            'separator2':['separator', None],
                                            'nonbonded': ['MenuItem', None],
                    
                                           }
                               ],
                    
                    
                    'hide'   : [
                                'submenu',  {
                                            'lines'    : ['MenuItem', menu_hide_lines],
                                            'sticks'   : ['MenuItem', menu_hide_sticks],
                                            'spheres'  : ['MenuItem', menu_hide_spheres],
                                            'nonbonded': ['MenuItem', None],
                                            }
                                ],
                    
                    
                    'separator2':['separator', None],

                    
                    
                    'label':  ['submenu' , {
                                            'Atom'         : [
                                                               'submenu', {
                                                                           'lines'    : ['MenuItem', None],
                                                                           'sticks'   : ['MenuItem', None],
                                                                           'spheres'  : ['MenuItem', None],
                                                                           'nonbonded': ['MenuItem', None],
                                                                           }
                                                              ],
                                            
                                            'atomic index' : ['MenuItem', None],
                                            'residue name' : ['MenuItem', None],
                                            'residue_index': ['MenuItem', None],
                                           },
                               ]
                    }



        if pick_menu is None:
            ''' Standard Sele Menu '''
            pick_menu = { 
                    'header' : ['MenuItem', None],
                    
                    
                    
                    'separator1':['separator', None],
                    
                    
                    'show'   : [
                                'submenu' ,{
                                            
                                            'lines'         : ['MenuItem', menu_show_lines],
                                            'sticks'        : ['MenuItem', menu_show_sticks],
                                            'spheres'       : ['MenuItem', menu_show_spheres],
                                            'dynamic bonds' : ['MenuItem', dynamic_test],
                                            'separator2'    : ['separator', None],
                                            'nonbonded'     : ['MenuItem', None],
                    
                                           }
                               ],
                    
                    
                    'hide'   : [
                                'submenu',  {
                                            'lines'    : ['MenuItem', menu_hide_lines],
                                            'sticks'   : ['MenuItem', menu_hide_sticks],
                                            'spheres'  : ['MenuItem', menu_hide_spheres],
                                            'nonbonded': ['MenuItem', None],
                                            }
                                ],
                    
                    
                    'separator2':['separator', None],

                    }




        self.glwidget.build_glmenu(bg_menu   = bg_menu, 
                                   sele_menu = sele_menu, 
                                   obj_menu  = obj_menu,
                                   pick_menu = pick_menu )

    



    def restart_session (self, filename = None):
        """ Function doc """
        
        self.main_session.window.destroy()
        self.main_session = None
        self.__init__(glwidget = True, toolkit = 'gtk3')
        self.treestore = Gtk.TreeStore(
                                                str , # Name
                                                
                                                bool, # toggle active=1
                                                bool, # radio  active=2
                                                
                                                bool, # toggle visible = 3
                                                bool, # radio  visible = 4
                                                
                                                int , # vismol_object index
                                                bool, # is vismol_object index visible?
                                                int , # pdynamo system index 
                                                bool, # is pdynamo system index visible?
                                                )
        self.parents = {}
        self.insert_glmenu()

        #self.main_session.window.destroy()
        window = EasyHybridMainWindow(self)


        #Gtk.main()

    def load_easyhybrid_data_to_session (self, easyhybrid_session_data):
        """ Function doc """
        
    
    def pDynamo_selections (self):
        """ Function doc """
        
        #atomref = AtomSelection.FromAtomPattern( self.cSystem, _centerAtom )
        #core    = AtomSelection.Within(self.cSystem,atomref,_radius)
        #core2   = AtomSelection.ByComponent(self.cSystem,core)
        #self.cSystem = PruneByAtom( self.cSystem , Selection(core2) )
        ##---------------------------------------------------
        #self.cSystem.label = self.baseName + "#{} Pruned System ".format(self.systemCoutCurr) 
        #self.cSystem.DefineNBModel( self.nbModel )
        #self.cSystem.Energy()





def main():
    
    #vismolSession  =  VisMolSession(glwidget = True, toolkit = 'gtk3')
    vismolSession  =  EasyHybridVismolSession(glwidget = True, toolkit = 'gtk3')
    
    vismolSession.treestore = Gtk.TreeStore(
                                            str  ,   #                                   # 0
                                            bool ,   # toggle active=1                   # 1
                                            bool ,   # toggle visible = 3                # 2 
                                                                                         
                                            bool ,   # radio  active  = 2                # 3      
                                            bool ,   # radio  visible = 4                # 4     
                                                                                         
                                            bool  ,  # traj radio  active = 5            # 5        
                                            bool  ,  # is trajectory radio visible?      # 6          
                                                                                         
                                            int,     #                                   # 7
                                            int,     # pdynamo system index              # 8
                                            int,)    # frames  # 9
    
    
    
    '''
    vismolSession.treestore = Gtk.TreeStore(
                                        str , # Name
                                        
                                        bool, # toggle active=1
                                        bool, # radio  active=2
                                        
                                        bool, # toggle visible = 3
                                        bool, # radio  visible = 4
                                        
                                        bool, # radio  active = 5 
                                        #int , # vismol_object index
                                        bool, # is trajectory toogle visible?
                                        int , # pdynamo system index 
                                        bool, # is pdynamo system index visible?
                                        )
    '''
    vismolSession.parents = {}
    vismolSession.insert_glmenu()
    window = EasyHybridMainWindow(vismolSession)
    Gtk.main()
 

if __name__ == '__main__':
    main()

