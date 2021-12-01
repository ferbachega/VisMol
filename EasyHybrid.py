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




class EasyHybridVismolSession(VisMolSession):
    """ Class doc """
    
    def add_vismol_object_to_vismol_session (self, pdynamo_session = None, 
                                                   rep             = True, 
                                                   vismol_object   = None, 
                                                   autocenter      = True):
        """ Function doc """
        vismol_object.index = self.vobj_counter
        self.vismol_objects.append(vismol_object)
        self.vismol_objects_dic[self.vobj_counter] = vismol_object
        self.vobj_counter += 1
        
        
        ##print('\n\n\n')
        #print (self.parents)
        #vobj_index = self.vismol_objects.index(vismol_object)
        vobj_index = vismol_object.index
        sys_index  = vismol_object.eb_sys_id
        
        
        if sys_index in self.parents.keys():
            pass
        else:          
            # Creates a new "parent" when a new system is loaded into memory. 
            for row in self.treestore:
                #row[2] = row.path == selected_path
                row[2] =  False
                for i,j in enumerate(row):
                    print(i, j,)
            
            self.parents[sys_index] = self.treestore.append(None, [pdynamo_session.systems[sys_index]['name'], False, True, False, True, vobj_index, vismol_object.eb_sys_id])

        self.treestore.append(self.parents[vismol_object.eb_sys_id], [vismol_object.name, True , False , True, False, vobj_index, vismol_object.eb_sys_id])

        print('\n\n\n')

        for vobj in self.vismol_objects:
            print(vobj.name, vobj.eb_sys_id)
        print('\n\n\n')
        print(vismol_object.index)

        if rep:
            #self.vismol_objects[-1].generate_indexesresentations (reps_list = self.indexes)
            #print (self.vismol_objects[-1].representations)

            #rep =  CartoonRepresentation(name = 'cartoon', active = True, _type = 'mol', visObj = self.vismol_objects[-1], glCore = self.glwidget.vm_widget)
            #self.vismol_objects[-1].representations[rep.name] = rep
            
            #rep =  RibbonsRepresentation(name = 'ribbons', active = True, _type = 'mol', visObj = self.vismol_objects[-1], glCore = self.glwidget.vm_widget)
            #self.vismol_objects[-1].representations[rep.name] = rep
            
            self.vismol_objects[-1].create_new_representation (rtype = 'lines')
            #rep  = LinesRepresentation (name = 'lines', active = True, _type = 'mol', visObj = self.vismol_objects[-1], glCore = self.glwidget.vm_widget)
            #self.vismol_objects[-1].representations[rep.name] = rep
            #
            rep  = NonBondedRepresentation (name = 'nonbonded', active = True, _type = 'mol', visObj = self.vismol_objects[-1], glCore = self.glwidget.vm_widget)
            self.vismol_objects[-1].representations[rep.name] = rep
        
            if autocenter:
                #print(self.vismol_objects[-1].mass_center)
                self.glwidget.vm_widget.center_on_coordinates(self.vismol_objects[-1], self.vismol_objects[-1].mass_center)
            else:
                self.glwidget.vm_widget.queue_draw()
            self.gtk_widgets_update ()
        






    
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
                selection = self.selections[self.current_selection]

                pdmsys_active = self.main_session.pDynamo_session.active_id
                self.main_session.pDynamo_session.systems[pdmsys_active]['qc_table'] = []
                
                if len(selection.selected_objects) >= 2:
                    print ('''The list of selected atoms contains items coming from different structures. These selections can present problems regarding the numbering of atoms (it is recommended that you select atoms from a single structure).
                    Do you want to continue anyway?''')
                
                elif selection.selected_atoms[0].Vobject.eb_sys_id != self.main_session.pDynamo_session.active_id:
                    print('The selected atoms do not belong to the active system in pdynamo.')
                
                
                for atom in selection.selected_atoms:
                    #print(atom.index-1, atom.name, atom.resn)
                    self.main_session.pDynamo_session.systems[pdmsys_active]['qc_table'].append(atom.index -1)
                #print('selection_qc',self.main_session.pDynamo_session.systems[pdmsys_active]['qc_table'] )
                self.main_session.run_dialog_set_QC_atoms()

            def set_as_free_atoms (_):
                """ Function doc """
                selection         = self.selections[self.current_selection]
                
                # these are the new atoms to bet set as fixed
                freelist = []                
                for atom in selection.selected_atoms:
                    #print(atom.index-1, atom.name, atom.resn)
                    freelist.append(atom.index -1)
                    atom.get_color()  
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
                    self.glwidget.vm_widget.queue_draw()
                #self.main_session.pDynamo_session.vismol_selection_qc = selection.copy()

            def set_as_fixed_atoms (_):
                """ Function doc """
                selection         = self.selections[self.current_selection]
                
                # these are the new atoms to bet set as fixed
                fixedlist = []                
                for atom in selection.selected_atoms:
                    #print(atom.index-1, atom.name, atom.resn)
                    fixedlist.append(atom.index -1)
                ##print('selection_free',fixedlist )
                #----------------------------------------------
                
                #Combining with list that the already exists
                pdmsys_active =   self.main_session.pDynamo_session.active_id
                
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

    












def main():
    
    #vismolSession  =  VisMolSession(glwidget = True, toolkit = 'gtk3')
    vismolSession  =  EasyHybridVismolSession(glwidget = True, toolkit = 'gtk3')
    vismolSession.treestore = Gtk.TreeStore(
                                            str , # Name
                                            
                                            bool, # toggle active=1
                                            bool, # radio  active=2
                                            
                                            bool, # toggle visible = 3
                                            bool, # radio  visible = 4
                                            
                                            int , # vismol_object index
                                            int , # pdynamo system index 
                                            
                                            
                                            )
                                            
   
    vismolSession.parents   = {}
    #vismolSession.vobj_counter  = 0
    
    
    
    
    '''
    def menu_show_lines (_):
        """ Function doc """
        vismolSession.show_or_hide( _type = 'lines', show = True)

    def menu_hide_lines (_):
        """ Function doc """
        vismolSession.show_or_hide( _type = 'lines', show = False)

    def menu_show_sticks (_):
        """ Function doc """
        vismolSession.show_or_hide( _type = 'sticks', show = True)

    def menu_hide_sticks (_):
        """ Function doc """
        vismolSession.show_or_hide( _type = 'sticks', show = False)

    def menu_show_spheres (_):
        """ Function doc """
        vismolSession.show_or_hide( _type = 'spheres', show = True)

    def menu_hide_spheres (_):
        """ Function doc """
        vismolSession.show_or_hide( _type = 'spheres', show = False)
    menu = { 
            'Teste Menu from main ' : ['MenuItem', None],
            
            
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

            }

    '''
    vismolSession.insert_glmenu()
    window = EasyHybridMainWindow(vismolSession)
    #window.connect("destroy", Gtk.main_quit)
    #window.show_all()
    Gtk.main()
    
    
    #print ('ops')
    ##vismolSession.insert_glmenu(bg_menu = menu)
    #vismolSession.insert_glmenu()
    #gui            = MainWindow(vismolSession)
    #gui.run()
    #print ('ops')

if __name__ == '__main__':
    main()

