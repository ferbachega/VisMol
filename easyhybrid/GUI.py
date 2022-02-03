import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
from GTKGUI.gtkWidgets.filechooser import FileChooser
from easyhybrid.pDynamoMethods.pDynamo2Vismol import *
import gc
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
        #print('self.Visible',self.Visible)
    
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
        #print(files)
    
    def on_import_files_button_clicked (self, button):
        """ Function doc """
        files = self.easyhybrid_main.filechooser.open(select_multiple = True)
        #print(files)

        for _file in files:
            #for res in self.VObj.chains[chain].residues:
                ##print(res.resi, res.resn, chain,  len(res.atoms) ) 
            systemtype = self.system_types_combo.get_active()
            filetype = self.filetype_parser( _file, systemtype)
            self.residue_liststore.append(list([_file, filetype, '10' ]))
        self.treeview.set_model(self.residue_liststore)
        self.files['opls_folder'] =  self.builder.get_object('OPLS_folderchooserbutton').get_filename()
        #print(self.files)


    def on_button_import_a_new_system_clicked (self, button):
        """ Function doc """
        
        if button == self.builder.get_object('ok_button_import_a_new_system'):
            print('ok_button_import_a_new_system')
            #self.on_button1_clicked_create_new_project(button)
            #self.dialog.hide()
        if button == self.builder.get_object('cancel_button_import_a_new_system'):
            print('cancel_button_import_a_new_system')
            self.dialog.hide()
            
    def on_button4_import_system_clicked (self, button):
        #print('ok_button_import_a_new_system')
        systemtype = self.system_types_combo.get_active()
        
        name =  self.builder.get_object('entry_system_name').get_text()

        self.easyhybrid_main.pDynamo_session.load_a_new_pDynamo_system_from_dict(filesin = self.files, 
                                                                                 systype = systemtype, 
                                                                                 name = name)
        
        #if systemtype == 2:
        #    self.files['opls_folder'] =  self.builder.get_object('OPLS_folderchooserbutton').get_filename()
        
        #print ('systemtype',systemtype, self.files )
        ##self.easyhybrid_main.pDynamo_session.get_bonds_from_pDynamo_system()
        #
        #vismol_object = self.easyhybrid_main.pDynamo_session.build_vismol_object_from_pDynamo_system (name = name)
        self.CloseWindow(button, data  = None)


class EasyHybridImportTrajectoryWindow:
    """ Class doc """
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/easyhybrid_import_trajectory_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('import_trajectory_window')
            self.window.set_title('Import Trajectory Window')
            self.window.set_keep_above(True)
            '''--------------------------------------------------------------------------------------------'''
            
            self.combobox_pdynamo_system = self.builder.get_object('combobox_pdynamo_system')
            renderer_text = Gtk.CellRendererText()
            self.combobox_pdynamo_system.add_attribute(renderer_text, "text", 0)

            self.system_liststore = Gtk.ListStore(str)
            
            names = [ ]
            for system_key in self.easyhybrid_main.pDynamo_session.systems.keys():
                sys = self.easyhybrid_main.pDynamo_session.systems[system_key]
                names.append(sys['name'] )
                print(sys['name'] )
            
            for name in names:
                self.system_liststore.append([name])
            
            self.combobox_pdynamo_system.set_model(self.system_liststore)
            self.combobox_pdynamo_system.set_active(0)
            
            #methods = [
            #    "Conjugate Gradient" ,
            #    "FIRE"               ,
            #    "L-BFGS"             ,
            #    "Steepest Descent"   ,
            #    ]
            #for method in methods:
            #    self.method_store.append([method])
            #    print (method)
            #
            #self.methods_combo = self.builder.get_object('combobox_geo_opt')
            #self.methods_combo.set_model(self.method_store)
            #self.methods_combo.connect("changed", self.on_name_combo_changed)
            #self.methods_combo.set_model(self.method_store)
            #
            #renderer_text = Gtk.CellRendererText()
            #self.methods_combo.pack_start(renderer_text, True)
            #self.methods_combo.add_attribute(renderer_text, "text", 0)
            #'''--------------------------------------------------------------------------------------------'''
            self.combox = self.builder.get_object('combobox_trajectory_type')
            self.combox.set_active(0)
            self.window.show_all()
            self.Visible  = True
    
    def import_data (self, button):
        """ Function doc """
        entry_name    = None
        idnum     = self.combobox_pdynamo_system.get_active()
        text      = self.combobox_pdynamo_system.get_active_text()
        
        print(idnum, text )
        #traj_log      = int  ( self.builder.get_object('entry_traj_log').get_text() )
        #logFrequency  = int  ( self.builder.get_object('entry_log_frequence').get_text())
        #max_int       = int  ( self.builder.get_object('entry_max_int').get_text()  )
        #rmsd_tol      = float( self.builder.get_object('entry_rmsd_tol').get_text() )
        #entry_name    =        self.builder.get_object('entry_name').get_text()  
        #
        #save_trajectory = self.builder.get_object('checkbox_save_traj').get_active() 
    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
    
    
    def __init__(self, main = None):
        """ Class initialiser """
        self.easyhybrid_main     = main
        self.Visible             =  False        




class EasyHybridGeometryOptimizatrionWindow(Gtk.Window):
    """ Class doc """
    
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/easyhybrid_geometry_optimization_window.glade'))
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
            '''--------------------------------------------------------------------------------------------'''
            self.methods_combo.set_active(0)
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


    def run_opt (self, button):
        """ Function doc """
        entry_name    = None
        method_id     = self.builder.get_object('combobox_geo_opt').get_active()
        traj_log      = int  ( self.builder.get_object('entry_traj_log').get_text() )
        logFrequency  = int  ( self.builder.get_object('entry_log_frequence').get_text())
        max_int       = int  ( self.builder.get_object('entry_max_int').get_text()  )
        rmsd_tol      = float( self.builder.get_object('entry_rmsd_tol').get_text() )
        entry_name    =        self.builder.get_object('entry_name').get_text()  
        
        save_trajectory = self.builder.get_object('checkbox_save_traj').get_active() 
        
        if method_id == 0:
            self.easyhybrid_main.pDynamo_session.run_ConjugateGradientMinimize_SystemGeometry(                  
                                                                                           
                                                                                           logFrequency         = logFrequency  , 
                                                                                           
                                                                                           maximumIterations    = max_int       , 
                                                                                           
                                                                                           rmsGradientTolerance = rmsd_tol      , 
                                                                                           save_trajectory      = save_trajectory,
                                                                                           trajectory_path      = None  
                                                                                           )


    
        self.window.destroy()
        self.Visible    =  False
    
    def on_name_combo_changed(self, widget):
        """ Function doc """
        print('eba - apagar')

    
    def _ (_):
        """ Function doc """
        
   
   
   

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
        #print('self.Visible',self.Visible)
    
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
        #print(button)
        #charge         = self.spinbutton_charge.get_value_as_int()
        #multiplicity   = self.spinbutton_multiplicity.get_value_as_int()
        #print('\n\ncharge'  , self.charge      )
        #print('multiplicity', self.multiplicity)
        #print('method_id'   , self.method_id   )
        
        if self.builder.get_object('radio_button_restricted').get_active():
            #print("%s is active" % (self.builder.get_object('radio_button_restricted').get_label()))
            self.restricted = True
        else:
            #print("%s is not active" % (self.builder.get_object('radio_button_restricted').get_label()))
            self.restricted = False
        
        
        parameters = {
                    'charge'       : self.charge      ,
                    'multiplicity' : self.multiplicity,
                    'method'       : self.methods_id_dictionary[self.method_id]   ,
                    'restricted'   : self.restricted  ,
                    
                     
                     }
        
        
        ##print(parameters)
        
        self.easyhybrid_main.pDynamo_session.define_a_new_QCModel(parameters =parameters)
        
        
        #self.easyhybrid_main.vismolSession.
        
        self.window.destroy()
        self.Visible    =  False


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



class EasyHybridMergeSystem(Gtk.Window):
    """ Class doc """
    
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/merge_systems.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('Merge_systems_window')
            self.box1   = self.builder.get_object('box1')
            self.box2   = self.builder.get_object('box2')
            self.entry  = self.builder.get_object('merge_entry')
            '''--------------------------------------------------------------------------------------------'''
            self.system_type_store = Gtk.ListStore(int, str)
            
            for index, system in self.easyhybrid_main.pDynamo_session.systems.items():
                
                name  = "{} {}".format(index, system['name'] )
                print(name)
                #name = 'teste'
                self.system_type_store.append([index, name])
           
            
            self.combo1 = Gtk.ComboBox.new_with_model(self.system_type_store)
            self.combo2 = Gtk.ComboBox.new_with_model(self.system_type_store)
            
            self.box1.pack_start(self.combo1, True, True, 0)
            self.box2.pack_start(self.combo2, True, True, 0)

            
            renderer_text = Gtk.CellRendererText()
            self.combo1.pack_start(renderer_text, True)
            self.combo1.add_attribute(renderer_text, "text", 1)
            
            renderer_text = Gtk.CellRendererText()
            self.combo2.pack_start(renderer_text, True)
            self.combo2.add_attribute(renderer_text, "text", 1)
            '''--------------------------------------------------------------------------------------------'''

            
            self.window.show_all()                                               
            self.builder.connect_signals(self)                                   

            self.Visible  =  True
            #----------------------------------------------------------------
    
    def ok_button (self, button):
        """ Function doc """
        print(button)
        print(self.combo1.get_active())
        print(self.combo1.get_active_id())
        print(self.combo1.get_active_iter())
        print(self.combo2.get_active())
        
        name1 = None
        name2 = None
        
        tree_iter = self.combo1.get_active_iter()
        if tree_iter is not None:
            model    = self.combo1.get_model()
            name1    = model[tree_iter][1]
            index1   = model[tree_iter][0]
            print("Selected: system =%s" % index1, name1)
        
        tree_iter = self.combo2.get_active_iter()
        if tree_iter is not None:
            model    = self.combo2.get_model()
            name2    = model[tree_iter][1]
            index2   = model[tree_iter][0]
            print("Selected: system =%s" % index2, name2)
        
        
        new_system_label = self.entry.get_text()
        
        
        if index2 != index1 and name1 is not None and name2 is not None:
            
            
            system1 = self.easyhybrid_main.pDynamo_session.systems[index1]['system']
            system2 = self.easyhybrid_main.pDynamo_session.systems[index2]['system']
            system1.Summary()
            system2.Summary()
            print(system1)
            print(system2)
            self.easyhybrid_main.pDynamo_session.merge_systems (system1 = system1, 
                                                                system2 = system2, 
                                                                label   = 'Merged System', 
                                                                summary = True)
            
            
            
            
            
        
        
        
    
    
    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
    
    def __init__(self, main = None):
        """ Class initialiser """
        self.easyhybrid_main     = main
        self.Visible             =  False        
        





#
#class EasyHybridMergeSystem:
#    """ Class doc """
#    def OpenWindow (self):
#        """ Function doc """
#        if self.Visible  ==  False:
#            self.builder = Gtk.Builder()
#            self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/merge_systems.glade'))
#            self.builder.connect_signals(self)
#            
#            self.window = self.builder.get_object('Merge_systems_window')
#            self.window.set_title('Merge pDynamo Systems - Window')
#            self.window.set_keep_above(True)
#            '''--------------------------------------------------------------------------------------------'''
#            
#            #'''
#            self.system_liststore = Gtk.ListStore(int, str)
#            
#            self.box1 = self.builder.get_object('box1')
#            self.box2 = self.builder.get_object('box2')
#            self.main_box = self.builder.get_object('main_box')
#            #renderer_text  = Gtk.CellRendererText()
#            #renderer_text2 = Gtk.CellRendererText()
#            
#            #self.combobox_1.add_attribute(renderer_text, "text", 1)
#            #self.combobox_2.add_attribute(renderer_text2, "text", 1)
#            
#            
#            for index, system in self.easyhybrid_main.pDynamo_session.systems.items():
#                
#                name  = "{} {}".format(index, system['name'] )
#                print(name)
#                name = 'teste'
#                self.system_liststore.append([1,name])
#
#            #self.system_types_combo = Gtk.ComboBox.new_with_model(self.system_type_store)
#            
#            
#            self.combobox_1 = Gtk.ComboBox.new_with_model(self.system_liststore) 
#            renderer_text = Gtk.CellRendererText()
#            self.combobox_1.pack_start(renderer_text, True)
#            
#            self.combobox_2 = Gtk.ComboBox.new_with_model(self.system_liststore)
#            renderer_text = Gtk.CellRendererText()            
#            self.combobox_2.pack_start(renderer_text, True)
#
#            
#            self.box1.pack_start(self.box1, True, True, 0)
#            self.box2.pack_start(self.box2, True, True, 0)
#            
#            #self.button = Gtk.Button()
#            #self.main_box.pack_start(self.button, True, True, 0)
#            
#            self.combobox_1.show_all()
#            self.combobox_2.show_all()
#
#            self.window.show_all()
#            self.Visible  = True
#    
#    def ok_button (self, button):
#        """ Function doc """
#        print(button)
#        #entry_name    = None
#        #idnum     = self.combobox_pdynamo_system.get_active()
#        #text      = self.combobox_pdynamo_system.get_active_text()
#        
#        #print(idnum, text )
#    
#    def CloseWindow (self, button, data  = None):
#        """ Function doc """
#        self.window.destroy()
#        self.Visible    =  False
#    
#    
#    def __init__(self, main = None):
#        """ Class initialiser """
#        self.easyhybrid_main     = main
#        self.Visible             =  False        
#


class EasyHybridMainWindow ( ):
    """ Class doc """

    def __init__ (self, vismolSession = None):
        """ Class initialiser """
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(VISMOL_HOME,'GTKGUI/MainWindow.glade'))
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
        #self.treeview = GtkMainTreeView(vismolSession)
        self.treeview = GtkEasyHybridMainTreeView(self, vismolSession)
        
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
        self.NewSystemWindow              = EasyHybridImportANewSystemWindow(main = self)
        self.setup_QCModel_window         = EasyHybridSetupQCModelWindow(main = self)
        self.import_trajectory_window     = EasyHybridImportTrajectoryWindow(main = self)
        self.geometry_optimization_window = EasyHybridGeometryOptimizatrionWindow(main = self)
        self.merge_pdynamo_systems_window = EasyHybridMergeSystem(main = self)
        
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
                self.vismolSession.load_easyhybrid_serialization_file(filename)
                #infile = open(filename,'wb')
                #data = pickle.load(infile)
                #self.vismolSession.
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
            self.vismolSession.save_serialization_file()

            
            
            

        if button  == self.builder.get_object('toolbutton_energy'):
            energy = self.pDynamo_session.get_energy()
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
        
        if button  == self.builder.get_object('run_md'):
            #self.geometry_optimization_window.OpenWindow()
            #self.pDynamo_session.run_ConjugateGradientMinimize_SystemGeometry()
            self.import_trajectory_window.OpenWindow()
            #self.pDynamo_session.import_trajectory()
    
    
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

        #print("was turned", state)            
    
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
        print(menuitem)
        
        if menuitem == self.builder.get_object('menu_item_merge_system'):
            print(menuitem, 'menu_item_merge_system')
            self.merge_pdynamo_systems_window.OpenWindow()

class GtkEasyHybridMainTreeView(Gtk.TreeView):
    
    def __init__ (self,  main, vismolSession):
        Gtk.TreeView.__init__(self)
        self.main_session  = main
        self.vismolSession = vismolSession
        self.treeview_menu = TreeViewMenu(self)
        self.treestore     = self.vismolSession.treestore
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


    def on_cell_toggled(self, widget, path):
        
        self.treestore[path][1] = not self.treestore[path][1]
        #for i in path:
        #print(self.treestore[path][1], path, self.treestore[path][0],self.treestore[path][-1] )


        if self.treestore[path][1]:
            obj_index = self.treestore[path][7]
            #print('obj_index', obj_index)
            self.vismolSession.enable_by_index(index = int(obj_index))
            self.vismolSession.glwidget.queue_draw()
        else:
            obj_index = self.treestore[path][7]

            self.vismolSession.disable_by_index(index = int(obj_index))
            self.vismolSession.glwidget.queue_draw()


    def on_cell_radio_toggled2(self, widget, path):
        #widget.set_active()

        #selected_path = Gtk.TreePath(path)
        #print('selected_path:', selected_path)
        
        print('path:', path)
        print(widget)
        
        for treeview_iter in self.vismolSession.gtk_treeview_iters:
            self.treestore[treeview_iter][5] = False
            print(self.treestore[treeview_iter][0])
        self.treestore[path][5] = True
        #paths = []
        
        #for i, row in enumerate(self.treestore):
        #    #print (i, row, row.path)
        #    self.treestore[row.path][5] = False
        #    paths.append(row.path)
        #    #row[5] = False #row.path == selected_path
        
        #print ( '\n path:              ', path, 
        #        '\n path:              ', type(path), 
        #        '\n treestore:         ', self.treestore['0:2'][5], 
        #        '\n treestore[path][5]:', self.treestore)
        
        
        '''
        for row in self.treestore:
            print(row, self.treestore['0:2'])
            for row2 in row:
                print(row2)
        '''
        
  
    def on_cell_radio_toggled(self, widget, path):
        selected_path = Gtk.TreePath(path)
        print('selected_path', selected_path)
        print(widget)
        
        for row in self.treestore:
            row[3] = row.path == selected_path
            if row[3]:
                self.main_session.pDynamo_session.active_id = row[8]
            else:
                pass
            
            #for i,j in enumerate(row):
            #    print(i, j, 'row[2]', row[2], row[8],selected_path,row.path)#(row[2], row.path, selected_path)
        
        print('\n\nactive_id', self.main_session.pDynamo_session.active_id,'\n\n')


    def on_treeview_Objects_button_release_event(self, tree, event):
        if event.button == 3:
            
            
            selection     = self.get_selection()
            #print(selection)
            model         = self.get_model()
            #print(model)
            #print(self.treestore)
            
            (model, iter) = selection.get_selected()
            for item in model:
                print (item[0], model[iter][0])
            print (model[iter][:], iter, model, tree )
            if iter != None:
                self.selectedID  = str(model.get_value(iter, 1))  # @+
                self.selectedObj = str(model.get_value(iter, 2))
                print(self.selectedID, self.selectedObj)
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
            print(self.selectedID, model.get_value(iter, 7))
            print (model[iter][:], iter)
            visObj = self.vismolSession.vismol_objects_dic[self.selectedID]
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
                    'rename'                : self.f1 ,
                    'info'                  : self.f1 ,
                    'load data into system' : self.f1 ,
                    'define color palette'  : self.f2 ,
                    'edit parameters'       : self.f2 ,
                    'export as...'          : self.f3 ,
                    'merge system with...'  : self.f3 ,
                    'delete'                : self.f3 ,
                    #'test'  : self.f1 ,
                    #'f1'    : self.f1 ,
                    #'f2'    : self.f2 ,
                    #'gordão': self.f3 ,
                    #'delete': self.f3 ,
                    }
        self.build_tree_view_menu(functions)



    def f1 (self, visObj = None ):
        """ Function doc """
        selection        = self.treeview.get_selection()
        (model, iter)    = selection.get_selected()
        self.selectedID  = int(model.get_value(iter, 1))  # @+
        
        visObj = self.treeview.vismolSession.vismol_objects[self.selectedID]
        
        infile = self.filechooser.open()
        
        self.treeview.vismolSession.load_xyz_coords_to_vismol_object(infile , visObj)
        
        print (infile)
        
        
        self.treeview.store .clear()
        #self.vismolSession.vismol_objects_dic.items()
        for index, vis_object in self.treeview.vismolSession.vismol_objects_dic.items():
            print ('\n\n',vis_object.name,'\n\n')
            data = [vis_object.active          , 
                    #str(self.treeview.vismolSession.vismol_objects.index(vis_object)),
                    str(index),
                    vis_object.name            , 
                    str(len(vis_object.atoms)) , 
                    str(len(vis_object.frames)),
                   ]
            model.append(data)

    def f2 (self, visObj = None):
        """ Function doc """
        #print('f2')
        #self._show_lines(visObj = self.vismol_objects[0], indices = [0,1,2,3,4] )
        self.treeview.vismolSession.go_to_atom_window.OpenWindow()

    def f3 (self, visObj = None):
        """ Function doc """
        
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()


        self.selectedID  = int(model.get_value(iter, 1))  # @+
        
        
        
        del self.treeview.vismolSession.vismol_objects_dic[self.selectedID]
        '''
        visObj = self.treeview.vismolSession.vismol_objects_dic.pop(self.selectedID)
        del visObj
        '''
        self.treeview.store.clear()
        #n = 0
        #i = 1
        
        #self.vismolSession.vismol_objects_dic.items()
        #for vis_object in self.treeview.vismolSession.vismol_objects:
        for vobj_index ,vis_object in self.treeview.vismolSession.vismol_objects_dic.items():
            print ('\n\n',vis_object.name,'\n\n')

            data = [vis_object.active          , 
                    str(vobj_index),
                    vis_object.name            , 
                    str(len(vis_object.atoms)) , 
                    str(len(vis_object.frames)),
                   ]
            model.append(data)
        self.treeview.vismolSession.glwidget.queue_draw()
            #i +=1
            #n = n + 1
        
        
        #self.treeview.vismolSession.center(visObj)

        
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







