import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
from GTKGUI.gtkWidgets.filechooser import FileChooser
from easyhybrid.pDynamoMethods.pDynamo2Vismol import *

from easyhybrid.gui.easyhybrid_pDynamo_selection            import  PDynamoSelectionWindow
from easyhybrid.gui.ImportNewSystem                         import  ImportANewSystemWindow 
from easyhybrid.gui.easyhybrid_import_trajectory_window     import  ImportTrajectoryWindow 
from easyhybrid.gui.easyhybrid_export_systems               import  ExportDataWindow 
from easyhybrid.gui.PES_scan_window                         import  PotentialEnergyScanWindow 
from easyhybrid.gui.molecular_dynamics_window               import  MolecularDynamicsWindow 
from easyhybrid.gui.umbrella_sampling_window                import  UmbrellaSamplingWindow 
from easyhybrid.gui.geometry_optimization_window            import  GeometryOptimizatrionWindow 
from easyhybrid.gui.QCSetup_window                          import  EasyHybridSetupQCModelWindow 
from easyhybrid.gui.merge_systems                           import  MergeSystemsWindow 

import gc
import os

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')



class EasyHybridDialogGeneric(Gtk.Dialog):
    def __init__(self, parent, message):
        super().__init__(title="New QC list", transient_for=parent, flags=0)
        self.add_buttons(
            Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_YES, Gtk.ResponseType.YES
        )
        
        label = Gtk.Label(label=message)
        box = self.get_content_area()
        box.add(label)
        self.show_all()

        
class EasyHybridDialogSetQCAtoms(Gtk.Dialog):
    def __init__(self, parent):
        super().__init__(title="New QC list", transient_for=parent, flags=0)
        self.add_buttons(
            Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_YES, Gtk.ResponseType.YES
        )

        self.set_default_size(150, 100)

        label = Gtk.Label(label="A new quantum region has been defined. Would you like to set up your QC parameters now?")

        box = self.get_content_area()
        box.add(label)
        self.show_all()


class EasyHybridDialogEnergy(Gtk.Dialog):
    def __init__(self, parent, energy = None):
        super().__init__(title="Energy Dialog", transient_for=parent, flags=0)
        self.add_buttons(
             Gtk.STOCK_OK, Gtk.ResponseType.OK
        )

        self.set_default_size(300, 100)
        
        label = Gtk.Label(label="Energy = {0:.5f} (KJ/mol)".format(energy))

        box = self.get_content_area()
        box.add(label)
        self.show_all()




class EasyHybridMainWindow ( ):
    """ Class doc """

    def __init__ (self, vm_session = None):
        """ Class initialiser """
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/MainWindow.glade'))
        #self.builder.add_from_file(os.path.join(VISMOL_HOME,'GTKGUI/toolbar_builder.glade'))
        self.builder.connect_signals(self)
        self.window = self.builder.get_object('window1')
        self.window.set_default_size(1000, 600)                          
        
        #self.toolbar_builder = self.builder.get_object('toolbar_builder') 
        #self.builder.get_object('box1').pack_start(self.toolbar_builder, True, True, 1)
        
        # Status Bar
        self.statusbar_main = self.builder.get_object('statusbar1')
        #self.statusbar_main.push(0,'wellcome to EasyHydrid')
        self.statusbar_main.push(1,'welcome to EasyHybrid 3.0')
        
        self.paned_V = self.builder.get_object('paned_V')
        #self.nootbook  =  self.builder.get_object('notebook2')
        #self.window = Gtk.Window(title="VisMol window")
        #self.main_box = Gtk.Box()
        self.vm_session = vm_session#( main_session = None)
        self.vm_session.main_session = self
        
        self.window.connect("key-press-event",   self.vm_session.glwidget.key_pressed)
        self.window.connect("key-release-event", self.vm_session.glwidget.key_released)
        
                
        
        self.menu_box = self.builder.get_object('toolbutton_selection_box')
        self.box2 = self.builder.get_object('box2')
        self.selection_box = self.vm_session.selection_box
        #self.box2.pack_start(self.selection_box, True, True, 0)
        self.menu_box.add(self.selection_box)
        #remove this combobox for vismol tools after
        #self.combobox1 = self.builder.get_object('combobox1')
        #self.combobox1.set_model(self.vm_session.Vismol_selection_modes_ListStore)
        #self.renderer_text = Gtk.CellRendererText()
        #self.combobox1.pack_start(self.renderer_text, True)
        #self.combobox1.add_attribute(self.renderer_text, "text", 0)        
        '''This gtk list is declared in the VismolGLWidget file 
           (it does not depend on the creation of Treeview)'''
        #self.Vismol_Objects_ListStore = self.vm_session.Vismol_Objects_ListStore
        
        #self.box2.pack_start(self.toolbar_builder, True, True, 1)
        
        #-------------------------------------------------------------------      
        #                         notebook_V1
        #-------------------------------------------------------------------
        #self.notebook_V1 = Gtk.Notebook()
        #print (self.notebook_V1.set_tab_pos(Gtk.PositionType.LEFT))
        #self.page1 = Gtk.Box()
        #self.page1.set_border_width(5)
        
        #self.text_view = Gtk.TextView()
        #self.text_view.set_editable(True)
        #self.page1.add( self.text_view)
        
        #self.page1.add(Gtk.Label('Here is the content of the first section.'))
        #self.notebook_V1.append_page(self.page1, Gtk.Label('Logs'))
        
        #-------------------------------------------------------------------      
        #                         notebook_H1
        #-------------------------------------------------------------------
        self.notebook_H1 = Gtk.Notebook()
        self.ScrolledWindow_notebook_H1 = Gtk.ScrolledWindow()
        
        #self.Tree_notebook_H1           = Gtk.TreeView()
        #self.treeview = GtkMainTreeView(vm_session)
        self.treeview = GtkEasyHybridMainTreeView(self, vm_session)
        
        #self.treeview  = self.gtkTreeViewObj.treeview
        self.ScrolledWindow_notebook_H1.add(self.treeview)
        self.notebook_H1.append_page(self.ScrolledWindow_notebook_H1, Gtk.Label('Objects'))
        


        # the label we use to show the selection
        self.label = Gtk.Label()
        self.label.set_text("")

        
        #-------------------------------------------------------------------
        #                         notebook_H2
        #-------------------------------------------------------------------
        self.notebook_H2 = Gtk.Notebook()
        #-------------------------------------------------------------------
        
        self.paned_H = Gtk.Paned(orientation = Gtk.Orientation.HORIZONTAL)
        self.button = Gtk.Button(label="Click Here")
        #-------------------------------------------------------------------
        self.vm_session = vm_session#( main_session = None)
        self.filechooser   = FileChooser()
        #-------------------------------------------------------------------
        
        
        self.container = Gtk.Box (orientation = Gtk.Orientation.VERTICAL)
        self.command_line_entry = Gtk.Entry()

        
        if self.vm_session is not None:
            #player

            self.container.pack_start(self.vm_session.glwidget, True, True, 0)
            
            self.traj_frame = self.vm_session.trajectory_frame
            #self.container.pack_start(self.traj_frame, False, False, 1)
            #self.container.pack_start(self.command_line_entry, False, False, 0)

            self.notebook_H2.append_page(self.container, Gtk.Label('view'))
            self.notebook_H2.append_page(Gtk.TextView(), Gtk.Label('logs'))
            
            
            #self.HBOX = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 6)
            self.HBOX = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 0)
            self.HBOX.pack_start(self.notebook_H1, True, True, 0)
            self.HBOX.pack_start(self.traj_frame, False, False, 1)

            #self.paned_H.add(self.notebook_H1)
            self.paned_H.add(self.HBOX)
            self.paned_H.add(self.notebook_H2)

            self.paned_V.add(self.paned_H)
            #self.paned_V.add(Gtk.TextView())
            
            #self.paned_V.add(self.traj_frame)
        

        #self.player_frame = self.vm_session.player_frame
        #self.player_frame.show_all()
        
        
        '''#- - - - - - - - - - - -  pDynamo - - - - - - - - - - - - - - -#'''
        
        self.pDynamo_session = pDynamoSession(vm_session = vm_session)

        '''#- - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -#'''
        self.NewSystemWindow              = ImportANewSystemWindow(main = self)
        self.setup_QCModel_window         = EasyHybridSetupQCModelWindow(main = self)
        self.import_trajectory_window     = ImportTrajectoryWindow(main = self)
        self.geometry_optimization_window = GeometryOptimizatrionWindow(main = self)
        self.molecular_dynamics_window    = MolecularDynamicsWindow(main = self)
        self.merge_pdynamo_systems_window = MergeSystemsWindow(main = self)
        self.pDynamo_selection_window     = PDynamoSelectionWindow(main = self)
        self.PES_scan_window              = PotentialEnergyScanWindow(main=  self)
        self.umbrella_sampling_window     = UmbrellaSamplingWindow (main=  self)
        self.export_data_window           = ExportDataWindow(main=  self)
        
        self.window.connect("destroy", Gtk.main_quit)
        self.window.connect("delete-event",    Gtk.main_quit)
        self.window.show_all()

        #Gtk.main()

    def run (self):
        """ Function doc """
        Gtk.main()
    
    def run_dialog_set_QC_atoms (self, _type = None):
        """ Function doc """
        dialog = EasyHybridDialogSetQCAtoms(self.window)
        response = dialog.run()

        if response == Gtk.ResponseType.YES:
            #print("The OK button was clicked")
            self.setup_QCModel_window.OpenWindow()
        elif response == Gtk.ResponseType.CANCEL:
            
            print("The Cancel button was clicked")

        dialog.destroy()
    
    def gtk_load_files (self, button):
        filename = self.filechooser.open()
        #print('aqui ohh')
        

        
        if filename:
            if filename[-4:] == 'easy':
                print('ehf file')            
                self.vm_session.load_easyhybrid_serialization_file(filename)
                #infile = open(filename,'wb')
                #data = pickle.load(infile)
                #self.vm_session.
                #pickle.dump(data, outfile)
                #outfile.close()
            
            else:
                files = {'coordinates': filename}
                systemtype = 3
                self.pDynamo_session.load_a_new_pDynamo_system_from_dict(files, systemtype)
        else:
            pass

    def on_main_toolbar_clicked (self, button):
        """ Function doc """
        if button  == self.builder.get_object('toolbutton_new_system'):
            #self.dialog_import_a_new_systen = EasyHybridImportANewSystemDialog(self.pDynamo_session, self)
            #self.dialog_import_a_new_systen.run()
            #self.dialog_import_a_new_systen.hide()
            self.NewSystemWindow.OpenWindow()
        
        
        
        
        
        if button  == self.builder.get_object('toolbutton_save'):
            self.vm_session.save_serialization_file()

            
            
            

        if button  == self.builder.get_object('toolbutton_energy'):
            energy = self.pDynamo_session.get_energy()
            #self.pDynamo_session.charge_summary()
            
            #print(energy)
            dialog = EasyHybridDialogEnergy(parent = self.window, energy = energy)
            response = dialog.run()
            dialog.destroy()
            
        if button  == self.builder.get_object('toolbutton_setup_QCModel'):
            #self.dialog_import_a_new_systen = EasyHybridImportANewSystemDialog(self.pDynamo_session, self)
            #self.dialog_import_a_new_systen.run()
            #self.dialog_import_a_new_systen.hide()
            #self.NewSystemWindow.OpenWindow()
            self.setup_QCModel_window.OpenWindow()
        if button  == self.builder.get_object('toolbutton_geometry_optimization'):
            self.geometry_optimization_window.OpenWindow()
        
        if button  == self.builder.get_object('toolbutton_pDynamo_selections'):
            #print('toolbutton_pDynamo_selections')
            #print('self.pDynamo_session.picking_selections_list', self.vm_session.picking_selections.picking_selections_list)
            self.pDynamo_selection_window.OpenWindow()
            
            '''
            atom1 = self.vm_session.picking_selections.picking_selections_list[0]
            #print (atom1.chain, atom1.resn, atom1.resi, atom1.name)
            print ("%s:%s.%s:%s" %(atom1.chain, atom1.resn, atom1.resi, atom1.name))
            
            _centerAtom ="%s:%s.%s:%s" %(atom1.chain, atom1.resn, atom1.resi, atom1.name)
            _radius =  10.0
            
            self.pDynamo_session.selections (_centerAtom, _radius)
            '''
        if button  == self.builder.get_object('toolbutton_pes_scan'):
            self.PES_scan_window.OpenWindow()
        
        if button  == self.builder.get_object('toolbutton_molecular_dynamics'):
            self.molecular_dynamics_window.OpenWindow()

        if button  == self.builder.get_object('toolbutton_umbrella_sampling'):
            #print('toolbutton_umbrella_sampling')
            self.umbrella_sampling_window.OpenWindow()

            
    def menubar_togglebutton1 (self, button):
        """ Function doc """
        if button.get_active():
            state = "on"
            self.vm_session._picking_selection_mode = True
            button.set_label('Picking')
            
            
        else:
            state = "off"
            self.vm_session._picking_selection_mode = False
            button.set_label('Viewing')

    
    def test (self, widget):
        """ Function doc """
        container = Gtk.Box (orientation = Gtk.Orientation.VERTICAL)
        container.pack_start(self.notebook_V1, True, True, 0)
        container.pack_start(self.command_line_entry, False, False, 0)
        self.paned_V.add(container)
        self.paned_V.show()
        self.window.show_all()


    def on_toolbutton_trajectory_tool (self, button):
        """ Function doc """
        print (button)
        
    def on_main_menu_activate (self, menuitem):
        """ Function doc """
        #print(menuitem)
        
        if menuitem == self.builder.get_object('menu_item_merge_system'):
            #print(menuitem, 'menu_item_merge_system')
            self.merge_pdynamo_systems_window.OpenWindow()


    def update_gui_widgets (self, update_folder = True, update_coords = True):
        """ Function doc """
        
        # should be a function . in the future!!!
        if update_coords:
            #print(self.vm_session.starting_coords_liststore)
            starting_coords = []
            #self.easyhybrid_main.pDynamo_session
            self.vm_session.starting_coords_liststore.clear()
            for key, visObj in self.vm_session.vismol_objects_dic.items():
                
                #print(visObj.name, visObj.easyhybrid_system_id, visObj.active)
                
                if visObj.easyhybrid_system_id == self.pDynamo_session.active_id:
                    starting_coords.append([visObj.name, key])
            
            for item in starting_coords:
                self.vm_session.starting_coords_liststore.append(list(item))
                #print (item)
                
            
        
        '''--------------------------------------------------------------------------------------------'''
        if update_folder:
            active_id = self.pDynamo_session.active_id
            working_folder = self.pDynamo_session.systems[active_id]['working_folder']
            
            if self.geometry_optimization_window.Visible:
                self.geometry_optimization_window.update_working_folder_chooser(folder = working_folder)
            
            if self.molecular_dynamics_window.Visible:
                self.molecular_dynamics_window.update_working_folder_chooser(folder = working_folder)
        




class GtkEasyHybridMainTreeView(Gtk.TreeView):
    
    def __init__ (self,  main, vm_session):
        Gtk.TreeView.__init__(self)
        self.main_session  = main
        self.vm_session    = vm_session
        self.treeview_menu = TreeViewMenu(self)
        self.treestore     = self.vm_session.treestore
        #--------------------------------------------------------------                                  
                                                                                                         
        self.set_model(self.treestore)                                                                   
                                                                                                         
                                                                                                         
        #------------------ r a d i o  ------------------                                                
        renderer_radio = Gtk.CellRendererToggle()                                                        
        renderer_radio.set_radio(True)                                                                   
        renderer_radio.connect("toggled", self.on_cell_radio_toggled)                                    
        column_radio = Gtk.TreeViewColumn("active", renderer_radio,    active=3, visible = 4)            
        self.append_column(column_radio)                                                                 
                                                                                                         
                                                                                                         
        
        #------------------  t e x t  ------------------
        renderer_text = Gtk.CellRendererText()
        column_text = Gtk.TreeViewColumn("Object", renderer_text, text=0)
        column_text.set_sizing(Gtk.TreeViewColumnSizing.AUTOSIZE)
        column_text.set_resizable(True)
        self.append_column(column_text)  
        

        
        #----------------- t o g g l e ------------------
        renderer_toggle = Gtk.CellRendererToggle()
        renderer_toggle.connect("toggled", self.on_cell_toggled)
        #column_toggle = Gtk.TreeViewColumn("Visible", renderer_toggle, active=1, visible = 3)
        column_toggle = Gtk.TreeViewColumn("V", renderer_toggle, active=1, visible = 2)
        self.append_column(column_toggle)


        #------------------ r a d i o ------------------

        renderer_radio2 = Gtk.CellRendererToggle()
        renderer_radio2.set_radio(True)
        renderer_radio2.connect("toggled", self.on_cell_radio_toggled2)
        column_radio2 = Gtk.TreeViewColumn("T", renderer_radio2, active = 5, visible = 6)
        self.append_column(column_radio2) 


        #'''
        renderer_text = Gtk.CellRendererText()
        column_text = Gtk.TreeViewColumn("F", renderer_text, text=9,  visible = 6)
        self.append_column(column_text)  
        #'''


        self.connect('button-release-event', self.on_treeview_Objects_button_release_event )
        
        
        self.selected_path = None

    def on_cell_toggled(self, widget, path):
        
        self.treestore[path][1] = not self.treestore[path][1]
        #for i in path:
        #print(self.treestore[path][1], path, self.treestore[path][0],self.treestore[path][-1] )


        if self.treestore[path][1]:
            obj_index = self.treestore[path][7]
            #print('obj_index', obj_index)
            self.vm_session.enable_by_index(index = int(obj_index))
            self.vm_session.glwidget.queue_draw()
        else:
            obj_index = self.treestore[path][7]

            self.vm_session.disable_by_index(index = int(obj_index))
            self.vm_session.glwidget.queue_draw()


    def on_cell_radio_toggled2(self, widget, path):
        '''
        Change the radio button relative to the 'T' 
        collumn of the vismol objects treeview 
        '''

        #'''
        #print('\n\n\path:', path)
        #print('\n\n\path:', 'AQUI')
        #print(widget)
        
        for treeview_iter in self.vm_session.gtk_treeview_iters:
            self.treestore[treeview_iter][5] = False
            #print(self.treestore[treeview_iter][0],self.treestore[treeview_iter][-1],self.treestore[treeview_iter][-2] )
        self.treestore[path][5] = True
        
        
        vob_id = self.treestore[path][-3]
        size = len(self.main_session.vm_session.vismol_objects_dic[vob_id].frames)
        self.treestore[path][-1] = size
        size = self.treestore[path][-1]
        print (path, type(path),self.treestore['0:1'][0] ,  self.treestore[path][0], self.treestore[path][-3],  self.treestore[path][-1])
        self.main_session.vm_session.TrajectoryFrame.change_range(upper = size)
        
        #print('\n\n\path:', path)
        #'''
        
  
    def on_cell_radio_toggled(self, widget, path):
        '''This function changes the status of radio buttons in the main tree. 
        
        It must also change the working directory and the list of coordinates 
        that will be available in each system.
        
        see:
        
        self.vm_session.combobox_initial_coordinates
        self.vm_session.filechooser_working_folder  
        
        '''
    
        #print('\n\n\path:', 'AQUI')
        #print('\n\n\path:', path)
        #
        #selected_path = path
        self.selected_path = path
        selected_path      = Gtk.TreePath(path)
        
        #print(widget, selected_path)
        
        for row in self.treestore:
            row[3] = row.path == selected_path
            if row[3]:
                self.main_session.pDynamo_session.active_id = row[8]
            else:
                pass
        
        #print('\n\nactive_id', self.main_session.pDynamo_session.active_id,'\n\n')
        
        self.main_session.update_gui_widgets()

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    def on_treeview_Objects_button_release_event(self, tree, event):
        if event.button == 3:
            
            
            selection     = self.get_selection()
            #print(selection)
            model         = self.get_model()
            #print(model)
            #print(self.treestore)
            
            (model, iter) = selection.get_selected()
            for item in model:
                pass
                #print (item[0], model[iter][0])
            #print (model[iter][:], iter, model, tree )
            if iter != None:
                self.selectedID  = str(model.get_value(iter, 1))  # @+
                self.selectedObj = str(model.get_value(iter, 2))
                #print(self.selectedID, self.selectedObj)
                self.treeview_menu.open_menu(self.selectedObj)
                
                #self.builder.get_object('TreeViewObjLabel').set_label('- ' +self.selectedObj+' -' )
            
                #widget = self.builder.get_object('treeview_menu')
                #widget.popup(None, None, None, None, event.button, event.time)
                #print ('button == 3')


        if event.button == 2:
            selection     = tree.get_selection()
            model         = tree.get_model()
            (model, iter) = selection.get_selected()
            #pymol_object = model.get_value(iter, 0)
            #self.refresh_gtk_main_self.treeView()
            print ('button == 2')
            
            self.selectedID  = int(model.get_value(iter, 7))  # @+
            #print(self.selectedID, model.get_value(iter, 7))
            #print (model[iter][:], iter)
            visObj = self.vm_session.vismol_objects_dic[self.selectedID]
            self.vm_session.center(visObj)

        if event.button == 1:
            print ('event.button == 1:')
        
        
class TreeViewMenu:
    """ Class doc """
    
    def __init__ (self, treeview):
        """ Class initialiser """
        pass
        self.treeview = treeview
        self.filechooser   = FileChooser()
        functions = {
                    'rename'                : self.f2 ,
                    'info'                  : self.f2 ,
                    'load data into system' : self.load_data_to_a_system ,
                    'define color palette'  : self.f2 ,
                    'edit parameters'       : self.f2 ,
                    'export as...'          : self.menu_export_data_window ,
                    'merge system with...'  : self.f3 ,
                    'delete'                : self.f3 ,
                    #'test'  : self.f1 ,
                    #'f1'    : self.f1 ,
                    #'f2'    : self.f2 ,
                    #'gordão': self.f3 ,
                    #'delete': self.f3 ,
                    }
        self.build_tree_view_menu(functions)



    def menu_export_data_window (self,visObj = None ):
        """ Function doc """
        #self.treeview.main_session.export_data_window.OpenWindow(sys_selected = model.get_value(iter, 8))
        self.treeview.main_session.export_data_window.OpenWindow()
    def load_data_to_a_system (self, visObj = None ):
        """ Function doc """
        selection        = self.treeview.get_selection()
        (model, iter)    = selection.get_selected()
        self.selectedID  = int(model.get_value(iter, 1))  # @+
        
        #print(selection, model, iter, 'selectedID',self.selectedID)
        visObj = self.treeview.vm_session.vismol_objects_dic[self.selectedID]
        selection     = self.treeview.get_selection()
        
        (model, iter) = selection.get_selected()
        for item in model:
            pass
            #print (item[0], model[iter][0])
        #print (model[iter][:], iter, model, tree )
        if iter != None:
            self.selectedID  = str(model.get_value(iter, 1))  # @+
            self.selectedObj = str(model.get_value(iter, 2))
            print(self.selectedID, self.selectedObj, model.get_value(iter, 0),  model.get_value(iter, 8))
        self.treeview.main_session.import_trajectory_window.OpenWindow(sys_selected = model.get_value(iter, 8))


    def f2 (self, visObj = None):
        """ Function doc """
        #print('f2')
        #self._show_lines(visObj = self.vismol_objects[0], indices = [0,1,2,3,4] )
        self.treeview.vm_session.go_to_atom_window.OpenWindow()

    def f3 (self, visObj = None):
        """ Function doc """
        
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()


        self.selectedID  = int(model.get_value(iter, 1))  # @+
        
        
        
        del self.treeview.vm_session.vismol_objects_dic[self.selectedID]
        '''
        visObj = self.treeview.vm_session.vismol_objects_dic.pop(self.selectedID)
        del visObj
        '''
        self.treeview.store.clear()
        #n = 0
        #i = 1
        
        #self.vm_session.vismol_objects_dic.items()
        #for vis_object in self.treeview.vm_session.vismol_objects:
        for vobj_index ,vis_object in self.treeview.vm_session.vismol_objects_dic.items():
            print ('\n\n',vis_object.name,'\n\n')

            data = [vis_object.active          , 
                    str(vobj_index),
                    vis_object.name            , 
                    str(len(vis_object.atoms)) , 
                    str(len(vis_object.frames)),
                   ]
            model.append(data)
        self.treeview.vm_session.glwidget.queue_draw()
            #i +=1
            #n = n + 1
        
        
        #self.treeview.vm_session.center(visObj)

        
        #print('f3')

    def build_tree_view_menu (self, menu_items = None):
        """ Function doc """
        self.tree_view_menu = Gtk.Menu()
        for label in menu_items:
            mitem = Gtk.MenuItem(label)
            mitem.connect('activate', menu_items[label])
            self.tree_view_menu.append(mitem)
            #mitem = Gtk.SeparatorMenuItem()
            #self.tree_view_menu.append(mitem)

        self.tree_view_menu.show_all()

    def open_menu (self, visObj = None):
        """ Function doc """
        print (visObj)
        
        #print('AQ?UIIIIIIIIIIII')
        self.tree_view_menu.popup(None, None, None, None, 0, 0)




def check_filetype(filein):
    """ Function doc """
    
    data =  open(filein)







