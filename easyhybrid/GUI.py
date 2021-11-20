import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
from GTKGUI.gtkWidgets.filechooser import FileChooser
from easyhybrid.pDynamoMethods.pDynamo2Vismol import *


import os
VISMOL_HOME = os.environ.get('VISMOL_HOME')

#from GTKGUI.gtkWidgets.main_treeview import GtkMainTreeView


class EasyHybridImportANewSystemWindow(Gtk.Window):
    """ Class doc """
    
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/easyhybrid_widgets.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('ImportNewSystemWindow')
            self.window.set_border_width(10)
            self.window.set_default_size(500, 370)  

            
            '''--------------------------------------------------------------------------------------------'''
            self.system_type_store = Gtk.ListStore(str)
            system_types = [
                "AMBER",
                "CHARMM",
                "OPLS",
                "pdynamo files (*.pkl, *.yaml)",
                "other (*.pdb, *.xyz, *.mol2)",
                ]
            for system_type in system_types:
                self.system_type_store.append([system_type])
                print (system_type)
            
            self.system_types_combo = Gtk.ComboBox.new_with_model(self.system_type_store)
            #self.system_types_combo = self.builder.get_object('system_type_combox_from_import_a_new_system')
            self.box_combo = self.builder.get_object('box')
            self.box_combo.pack_start(self.system_types_combo, True, True, 0)
            
            self.system_types_combo.connect("changed", self.on_name_combo_changed)
            self.system_types_combo.set_model(self.system_type_store)
            
            renderer_text = Gtk.CellRendererText()
            self.system_types_combo.pack_start(renderer_text, True)
            self.system_types_combo.add_attribute(renderer_text, "text", 0)
            '''--------------------------------------------------------------------------------------------'''
           
            
            self.treeview = self.builder.get_object('gtktreeview_import_system')
            for i, column_title in enumerate(['file', "type", "number of atoms"]):
                renderer = Gtk.CellRendererText()
                column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                self.treeview.append_column(column)

            
            self.window.show_all()                                               
            self.builder.connect_signals(self)                                   
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()

            self.Visible  =  True
            
            self.files    = {
                            'amber_prmtop': None,
                            'charmm_par'  : [],
                            'charmm_psf'  : None,
                            'charmm_extra': None, 
                            'opls_folder' : [],
                            'coordinates' : None,
                            }
            self.system_types_combo.set_active(0)

            #----------------------------------------------------------------

    def CloseWindow (self, button, data  = None):
        """ Function doc """
        #self.BackUpWindowData()
        self.window.destroy()
        self.Visible    =  False
        print('self.Visible',self.Visible)
    
    def __init__(self, main = None):
        """ Class initialiser """
        self.easyhybrid_main     = main
        self.Visible             =  False        
        
        self.residue_liststore = Gtk.ListStore(str, str, str)
        #self.atom_liststore    = Gtk.ListStore(bool, int, str, str, int, int)
        #self.residue_filter    = False
        
        self.charmm_txt = '''For systems prepared in the traditional CHARMM / PSF format, the required files are: parameters (prm / par), topologies (psf) and coordinates (chm, crd, pdb, xyz, ...)

NOTE: You can include more than one parameter file if needed.'''
        self.amber_txt = '''Systems prepared using the AMBER force field require  two files: topologies (top / prmtop) and coordinates (crd, pdb, xyz, ...).  '''
        self.OPLS_txt = '''Systems prepared natively in pDynamo using the OPLS force field require: A parameter folder and a topology/coordinate file (pdb, mol2, mol, ...) .  '''
        self.gmx_txt = '''Systems prepared natively in GROMACS using the CHARMM(or AMBER) force field require: A parameter/topology (top) file and coordinate file (pdb, mol2, mol, ...) .  '''
        
    
    def on_name_combo_changed(self, widget):
        """ Function doc """
        fftype = self.system_types_combo.get_active()
        print (fftype)
        self.files    = {
            'amber_prmtop': None,
            'charmm_par'  : [],
            'charmm_psf'  : None,
            'charmm_extra': None, 
            'opls_folder' : [],
            'coordinates' : None,
            }
        self.residue_liststore = Gtk.ListStore(str, str, str)
        self.treeview.set_model(self.residue_liststore)
            
        if fftype == 0: #AMBER
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            self.builder.get_object('gtk_label_fftype').set_text(self.amber_txt)
 
            
        if fftype == 1: # "CHARMM":
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            self.builder.get_object('gtk_label_fftype').set_text(self.charmm_txt)

            
        if fftype == 10: #"GROMACS":
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            self.builder.get_object('gtk_label_fftype').set_text(self.gmx_txt)

            
        if fftype == 2:#"OPLS":
            self.builder.get_object('gtkbox_OPLS_folderchooser').show()
            self.builder.get_object('gtk_label_fftype').set_text(self.OPLS_txt)
            
        if fftype == 3: #"pDynamo files(*.pkl,*.yaml)":
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            
        if fftype == 4: #"Other(*.pdb,*.xyz,*.mol2...)":
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()


    def filetype_parser(self, filein, systemtype):
        filetype = self.GetFileType(filein)
        print (filetype, systemtype)
        if filetype in ['top', 'prmtop', 'TOP', 'PRMTOP', ]:
            if systemtype == 0:
                self.files['amber_prmtop'] = filein
                return 'amber parameters/topologies'
                
   
        
        elif filetype in ['par', 'prmtop', 'prm', 'PAR', 'PRM']:
            if systemtype == 1:
                self.files['charmm_par'].append(filein)
                return 'charmm parameters'
        
        
        
        
        elif filetype in ['psf', 'PSF','psfx', 'PSFX']:
            if systemtype == 1:
                self.files['charmm_psf'] = filein
                return 'charmm topologies'
        
        
        
        elif filetype in ['pdb', 'PDB','mol','MOL','mol2','MOL2', 'xyz', 'XYZ', 'crd', 'inpcrd', 'chm']:
            #if systemtype == 1:
            self.files['coordinates'] = filein
            return 'coordinates'
        
        
        
        elif filetype in ['pkl', 'PKL']:
            #if systemtype == 1:
            self.files['coordinates'] = filein
            return 'pDynamo coordinates'
        else:
            return 'unknow type'
        

    
    
    def GetFileType(self, filename):
        file_type = filename.split('.')
        return file_type[-1]



    def on_delete_files_button_clicked (self, button):
        """ Function doc """
        files = self.easyhybrid_main.filechooser.open(select_multiple = True)
        print(files)
    
    def on_import_files_button_clicked (self, button):
        """ Function doc """
        files = self.easyhybrid_main.filechooser.open(select_multiple = True)
        print(files)

        for _file in files:
            #for res in self.VObj.chains[chain].residues:
                #print(res.resi, res.resn, chain,  len(res.atoms) ) 
            systemtype = self.system_types_combo.get_active()
            filetype = self.filetype_parser( _file, systemtype)
            self.residue_liststore.append(list([_file, filetype, '10' ]))
        self.treeview.set_model(self.residue_liststore)
        self.files['opls_folder'] =  self.builder.get_object('OPLS_folderchooserbutton').get_filename()
        print(self.files)


    def on_button_import_a_new_system_clicked (self, button):
        """ Function doc """
        
        if button == self.builder.get_object('ok_button_import_a_new_system'):
            print('ok_button_import_a_new_system')
            #self.on_button1_clicked_create_new_project(button)
            #self.dialog.hide()
        if button == self.builder.get_object('cancel_button_import_a_new_system'):
            print('cancel_button_import_a_new_system')
            #self.dialog.hide()
    
    def on_button4_import_system_clicked (self, button):
        print('ok_button_import_a_new_system')
        systemtype = self.system_types_combo.get_active()
        self.easyhybrid_main.pDynamo_session.load_a_new_pDynamo_system_from_dict(self.files, systemtype)
        
        #if systemtype == 2:
        #    self.files['opls_folder'] =  self.builder.get_object('OPLS_folderchooserbutton').get_filename()
        
        print ('systemtype',systemtype, self.files )
        #self.easyhybrid_main.pDynamo_session.get_bonds_from_pDynamo_system()
        name =  self.builder.get_object('entry_system_name').get_text()
        
        vismol_object = self.easyhybrid_main.pDynamo_session.build_vismol_object_from_pDynamo_system (name = name)
        self.CloseWindow(button, data  = None)





class EasyHybridSetupQCModelWindow:
    """ Class doc """
    
    def __init__(self, main = None):
        """ Class initialiser """
        self.easyhybrid_main     = main
        self.Visible             =  False        
        
        self.methods_liststore = Gtk.ListStore(str, str, str)
        
        self.method_id    = 0    # 0 for am1, 1 for pm3  and...
        self.charge       = 0
        self.multiplicity = 1
        self.restricted   = True
        
        self.adjustment_charge = Gtk.Adjustment(value=0,
                                           lower=-100,
                                           upper=100,
                                           step_increment=1,
                                           page_increment=1,
                                           page_size=0)
        
        self.adjustment_multiplicity = Gtk.Adjustment(value=1,
                                           lower=1,
                                           upper=100,
                                           step_increment=1,
                                           page_increment=1,
                                           page_size=0)
        
        self.methods_id_dictionary = {
                                      0 : 'am1'     ,
                                      1 : 'am1dphot',
                                      2 : 'pm3'     ,
                                      3 : 'pm6'     ,
                                      4 : 'mndo'    ,
                                
                                      }
        
        
        
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/easyhybrid_QCSetup_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('SetupQCModelWindow')
            self.window.set_keep_above(True)
            #self.window.set_border_width(10)
            #self.window.set_default_size(500, 370)  

            
            '''--------------------------------------------------------------------------------------------'''
            self.methods_type_store = Gtk.ListStore(str)
            methods_types = [
                "am1",
                "am1dphot",
                "pm3",
                "pm6",
                "mndo",
                "ab initio - ORCA",
                "DFT / DFTB",
                ]
            for method_type in methods_types:
                self.methods_type_store.append([method_type])
                print (method_type)
            
            self.methods_combo = self.builder.get_object('QCModel_methods_combobox')
            self.methods_combo.connect("changed", self.on_name_combo_changed)
            self.methods_combo.set_model(self.methods_type_store)
            renderer_text = Gtk.CellRendererText()
            self.methods_combo.pack_start(renderer_text, True)
            self.methods_combo.add_attribute(renderer_text, "text", 0)
            self.methods_combo.set_active(self.method_id)
            '''--------------------------------------------------------------------------------------------'''


            self.spinbutton_charge       = self.builder.get_object('spinbutton_charge'      )
            self.spinbutton_multiplicity = self.builder.get_object('spinbutton_multiplicity')
            self.spinbutton_charge      .set_adjustment(self.adjustment_charge)
            self.spinbutton_multiplicity.set_adjustment(self.adjustment_multiplicity)
            
            self.window.show_all()                                               
            self.builder.connect_signals(self)                                   

            self.Visible  =  True
            #self.builder.get_object('dftp_setup_box').hide()
            #self.builder.get_object('orca_setup_box').hide()

    def CloseWindow (self, button, data  = None):
        """ Function doc """
        #self.BackUpWindowData()
        self.window.destroy()
        self.Visible    =  False
        print('self.Visible',self.Visible)
    
            #----------------------------------------------------------------
    def on_spian_button_change (self, widget):
        """ Function doc """
        self.charge       = self.spinbutton_charge.get_value_as_int()
        self.multiplicity = self.spinbutton_multiplicity.get_value_as_int()
        
        
    def on_name_combo_changed (self, combobox):
        """ Function doc """
        self.method_id = self.builder.get_object('QCModel_methods_combobox').get_active()
    
    def on_button_ok (self, button):
        """ Function doc """
        print(button)
        #charge         = self.spinbutton_charge.get_value_as_int()
        #multiplicity   = self.spinbutton_multiplicity.get_value_as_int()
        print('\n\ncharge'  , self.charge      )
        print('multiplicity', self.multiplicity)
        print('method_id'   , self.method_id   )
        
        if self.builder.get_object('radio_button_restricted').get_active():
            print("%s is active" % (self.builder.get_object('radio_button_restricted').get_label()))
            self.restricted = True
        else:
            print("%s is not active" % (self.builder.get_object('radio_button_restricted').get_label()))
            self.restricted = False
        
        
        parameters = {
                    'charge'       : self.charge      ,
                    'multiplicity' : self.multiplicity,
                    'method'       : self.methods_id_dictionary[self.method_id]   ,
                    'restricted'   : self.restricted  ,
                    
                     
                     }
        
        
        #print(parameters)
        
        self.easyhybrid_main.pDynamo_session.define_a_new_qcmodel(parameters =parameters)
        
        
        #self.easyhybrid_main.vismolSession.
        
        self.window.destroy()
        self.Visible    =  False
   
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

    def run_dialog_set_QC_atoms (self, _type = None):
        """ Function doc """
        dialog = EasyHybridDialogSetQCAtoms(self.window)
        response = dialog.run()

        if response == Gtk.ResponseType.YES:
            print("The OK button was clicked")
            self.setup_QCModel_window.OpenWindow()
        elif response == Gtk.ResponseType.CANCEL:
            print("The Cancel button was clicked")

        dialog.destroy()
    
    def gtk_load_files (self, button):
        filename = self.filechooser.open()
        
        if filename:
            files = {'coordinates': filename}
            systemtype = 3
            self.pDynamo_session.load_a_new_pDynamo_system_from_dict(files, systemtype)
           
            print ('systemtype',systemtype,files )
            name =  self.pDynamo_session.system.label
            
            vismol_object = self.pDynamo_session.build_vismol_object_from_pDynamo_system (name = name)
            
            
            
            
            #self.pDynamo_session.load_a_new_pDynamo_system (filein = filename)
            #self.pDynamo_session.get_bonds_from_pDynamo_system()
            #
            #vismol_object = self.pDynamo_session.build_vismol_object_from_pDynamo_system (name = filename)
            #vismol_object.set_model_matrix(self.vismolSession.glwidget.vm_widget.model_mat)
            #vismol_object.active = True
            #vismol_object._get_center_of_mass(frame = 0)
            #self.vismolSession.add_vismol_object_to_vismol_session (rep = True, vismol_object = vismol_object, autocenter =  True)
            
            
            #self.vismolSession.load(filename)
            
            #self.treeview.append()
            #self.treeview.refresh_gtk_main_treeview()
            #visObj = self.vismolSession.vismol_objects[-1]
            #self.treeview.append(visObj)

            #self.vismolSession.glwidget.vm_widget.center_on_coordinates(visObj, visObj.mass_center)
        else:
            pass


    def __init__ (self, vismolSession = None):
        """ Class initialiser """
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(VISMOL_HOME,'GTKGUI/MainWindow.glade'))
        self.builder.connect_signals(self)
        self.window = self.builder.get_object('window1')
        self.window.set_default_size(1000, 600)                          
        
        

        
        # togglebutton 
        #self.togglebutton1 = self.builder.get_object('togglebutton1')
        #self.togglebutton1.connect('clicked', self.menubar_togglebutton1)
        
        # Status Bar
        self.statusbar_main = self.builder.get_object('statusbar1')
        #self.statusbar_main.push(0,'wellcome to EasyHydrid')
        self.statusbar_main.push(1,'welcome to EasyHybrid 3.0')
        
        self.paned_V = self.builder.get_object('paned_V')
        #self.nootbook  =  self.builder.get_object('notebook2')
        #self.window = Gtk.Window(title="VisMol window")
        #self.main_box = Gtk.Box()
        self.vismolSession = vismolSession#( main_session = None)
        self.vismolSession.main_session = self
        
        self.window.connect("key-press-event",   self.vismolSession.glwidget.key_pressed)
        self.window.connect("key-release-event", self.vismolSession.glwidget.key_released)
        
                
        
        self.menu_box = self.builder.get_object('toolbutton_selection_box')
        self.box2 = self.builder.get_object('box2')
        self.selection_box = self.vismolSession.selection_box
        #self.box2.pack_start(self.selection_box, True, True, 0)
        self.menu_box.add(self.selection_box)
        #remove this combobox for vismol tools after
        #self.combobox1 = self.builder.get_object('combobox1')
        #self.combobox1.set_model(self.vismolSession.Vismol_selection_modes_ListStore)
        #self.renderer_text = Gtk.CellRendererText()
        #self.combobox1.pack_start(self.renderer_text, True)
        #self.combobox1.add_attribute(self.renderer_text, "text", 0)        
        '''This gtk list is declared in the VismolGLWidget file 
           (it does not depend on the creation of Treeview)'''
        #self.Vismol_Objects_ListStore = self.vismolSession.Vismol_Objects_ListStore
        
        
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
        self.treeview = GtkMainTreeView(vismolSession)
        
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
        self.vismolSession = vismolSession#( main_session = None)
        self.filechooser   = FileChooser()
        #-------------------------------------------------------------------
        
        
        self.container = Gtk.Box (orientation = Gtk.Orientation.VERTICAL)
        self.command_line_entry = Gtk.Entry()

        
        if self.vismolSession is not None:
            #player

            self.container.pack_start(self.vismolSession.glwidget, True, True, 0)
            
            self.traj_frame = self.vismolSession.trajectory_frame
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
        

        #self.player_frame = self.vismolSession.player_frame
        #self.player_frame.show_all()
        
        
        '''#- - - - - - - - - - - -  pDynamo - - - - - - - - - - - - - - -#'''
        
        self.pDynamo_session = pDynamoSession(vismolSession = vismolSession)

        '''#- - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -#'''
        self.NewSystemWindow      = EasyHybridImportANewSystemWindow(main = self)
        self.setup_QCModel_window = EasyHybridSetupQCModelWindow(main = self)

        self.window.connect("delete-event",    Gtk.main_quit)
        self.window.show_all()

        Gtk.main()

    def run (self):
        """ Function doc """
        Gtk.main()


    def on_main_toolbar_clicked (self, button):
        """ Function doc """
        if button  == self.builder.get_object('toolbutton_new_system'):
            #self.dialog_import_a_new_systen = EasyHybridImportANewSystemDialog(self.pDynamo_session, self)
            #self.dialog_import_a_new_systen.run()
            #self.dialog_import_a_new_systen.hide()
            self.NewSystemWindow.OpenWindow()
        
        if button  == self.builder.get_object('toolbutton_energy'):
            energy = self.pDynamo_session.get_energy()
            print(energy)
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
            self.pDynamo_session.run_ConjugateGradientMinimize_SystemGeometry()
            
            
    def menubar_togglebutton1 (self, button):
        """ Function doc """
        if button.get_active():
            state = "on"
            self.vismolSession._picking_selection_mode = True
            button.set_label('Picking')
            
            
        else:
            state = "off"
            self.vismolSession._picking_selection_mode = False
            button.set_label('Viewing')

        print("was turned", state)            
    
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
        




class GtkMainTreeView(Gtk.TreeView):
    """ Class doc """
    
    def __init__ (self, vismolSession):
        """ Class initialiser """
        
        Gtk.TreeView.__init__(self)
        self.vismolSession = vismolSession
        self.treeview_menu = TreeViewMenu(self)
        #self.store         = Gtk.ListStore(bool,str , str ,str, str)
        self.store         = vismolSession.Vismol_Objects_ListStore

        self.set_model(self.store)



        #----------------------------------------------------------------------
        # the cellrenderer for the second column - boolean rendered as a toggle
        renderer_toggle = Gtk.CellRendererToggle()
        # the second column is created
        column_in_out = Gtk.TreeViewColumn("", renderer_toggle, active=0)
        # and it is appended to the treeview
        self.append_column(column_in_out)
        # connect the cellrenderertoggle with a callback function
        renderer_toggle.connect("toggled", self.on_toggled)


        # the cellrenderer for text columns
        renderer_text = Gtk.CellRendererText()
        column = Gtk.TreeViewColumn("id",     renderer_text, text=1)
        self.append_column(column)

        renderer_text = Gtk.CellRendererText()
        column = Gtk.TreeViewColumn("Object", renderer_text, text=2)
        self.append_column(column)
    
        renderer_text = Gtk.CellRendererText()
        column = Gtk.TreeViewColumn("Atoms",  renderer_text, text=3)
        self.append_column(column)
    
        renderer_text = Gtk.CellRendererText()
        column = Gtk.TreeViewColumn("Frames", renderer_text, text=4)
        self.append_column(column)

        self.connect('button-release-event', self.on_treeview_Objects_button_release_event )
        #----------------------------------------------------------------------

    
    
    def on_toggled(self, widget, path):
        # the boolean value of the selected row
        current_value = self.store[path][0]

        # change the boolean value of the selected row in the model
        self.store[path][0] = not current_value       
        if self.store[path][0]:
            obj_index = self.store[path][1]
            self.vismolSession.enable_by_index(int(obj_index))
            self.vismolSession.glwidget.queue_draw()
        else:
            obj_index = self.store[path][1]
            self.vismolSession.disable_by_index(int(obj_index))
            self.vismolSession.glwidget.queue_draw()
 
    
    def append(self, visObj):
        i = self.vismolSession.vismol_objects.index(visObj)
        
        data = [visObj.active         , 
               str(i)                 ,
               visObj.name            , 
               str(len(visObj.atoms)) , 
               str(len(visObj.frames)),
               ]
        print (data)
        self.store.append(data)
        #self.set_model(self.liststore)
    
    def remove (self):
        """ Function doc """
        
    def on_treeview_Objects_button_release_event(self, tree, event):
        if event.button == 3:
            selection     = self.get_selection()
            model         = self.get_model()
            (model, iter) = selection.get_selected()
            if iter != None:
                self.selectedID  = str(model.get_value(iter, 1))  # @+
                self.selectedObj = str(model.get_value(iter, 2))
    
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
            
            self.selectedID  = int(model.get_value(iter, 1))  # @+
            visObj = self.vismolSession.vismol_objects[self.selectedID]
            self.vismolSession.center(visObj)

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
                'test':self.f1 ,
                'f1': self.f1,
                'f2': self.f2,
                'delete': self.f3,
        }
        self.build_glmenu(functions)



    def f1 (self, visObj = None ):
        """ Function doc """
        selection        = self.treeview.get_selection()
        (model, iter)    = selection.get_selected()
        self.selectedID  = int(model.get_value(iter, 1))  # @+
        
        visObj = self.treeview.vismolSession.vismol_objects[self.selectedID]
        
        infile = self.filechooser.open()
        
        self.treeview.vismolSession.load_xyz_coords_to_vismol_obejct(infile , visObj)
        
        print (infile)
        
        
        self.treeview.store .clear()
        for vis_object in self.treeview.vismolSession.vismol_objects:
            print ('\n\n',vis_object.name,'\n\n')
            data = [vis_object.active          , 
                    str(self.treeview.vismolSession.vismol_objects.index(vis_object)),
                    vis_object.name            , 
                    str(len(vis_object.atoms)) , 
                    str(len(vis_object.frames)),
                   ]
            model.append(data)
        #self.treeview.vismolSession.glwidget.queue_draw()
    
        #self.treeview.vismolSession.go_to_atom_window.OpenWindow()
    
    
    def f2 (self, visObj = None):
        """ Function doc """
        print('f2')
        #self._show_lines(visObj = self.vismol_objects[0], indices = [0,1,2,3,4] )
        self.treeview.vismolSession.go_to_atom_window.OpenWindow()

    def f3 (self, visObj = None):
        """ Function doc """
        
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()


        self.selectedID  = int(model.get_value(iter, 1))  # @+
        
        
        
        #visObj = self.treeview.vismolSession.vismol_objects[self.selectedID]
        visObj = self.treeview.vismolSession.vismol_objects.pop(self.selectedID)
        del visObj
        self.treeview.store .clear()
        #n = 0
        #i = 1
        for vis_object in self.treeview.vismolSession.vismol_objects:
            print ('\n\n',vis_object.name,'\n\n')

            data = [vis_object.active          , 
                    str(self.treeview.vismolSession.vismol_objects.index(vis_object)),
                    vis_object.name            , 
                    str(len(vis_object.atoms)) , 
                    str(len(vis_object.frames)),
                   ]
            model.append(data)
        self.treeview.vismolSession.glwidget.queue_draw()
            #i +=1
            #n = n + 1
        
        
        #self.treeview.vismolSession.center(visObj)

        
        print('f3')

    def build_glmenu (self, menu_items = None):
        """ Function doc """
        self.glMenu = Gtk.Menu()
        for label in menu_items:
            mitem = Gtk.MenuItem(label)
            mitem.connect('activate', menu_items[label])
            self.glMenu.append(mitem)
            mitem = Gtk.SeparatorMenuItem()
            self.glMenu.append(mitem)

        self.glMenu.show_all()

    def open_menu (self, visObj = None):
        """ Function doc """
        print (visObj)
        
        
        self.glMenu.popup(None, None, None, None, 0, 0)




def check_filetype(filein):
    """ Function doc """
    
    data =  open(filein)







