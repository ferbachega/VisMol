import pickle
from vModel import VismolObject
from OpenGL import GL
import ctypes

class LoadAndSaveFiles:
    """ Class doc """
    
    def __init__ (self, vm_session, pDynamo_session):
        """ Class initialiser """
        
        self.pDynamo_session = self.main_session.pDynamo_session
        self.main_session    = self.main_session
        
    def save_session (self, filename):
        """ Function doc """
  
        self.pDynamo_session = self.main_session.pDynamo_session
        
        vismol_objects_dic = {}
        
        #********************************************************************** 
        
        #---------------------------------------------------------------
        #                 V I S M O L    O B J E C T S
        #---------------------------------------------------------------
        
        for vobj_id, vobj in self.vismol_objects_dic.items():
            
            vismol_objects_dic[vobj_id] = {}
            '''-------------------------------------------------------'''
            '''Exporting atom list'''
            '''-------------------------------------------------------'''
            atoms = []
            for atom in vobj.atoms:
                
                atoms.append({
                              'index'      : atom.index-1    , 
                              'name'       : atom.name       , 
                              'resi'       : atom.resi       , 
                              'resn'       : atom.resn       , 
                              'chain'      : atom.chain      , 
                              'symbol'     : atom.symbol     , 
                              'occupancy'  : atom.occupancy  , 
                              'bfactor'    : atom.bfactor    , 
                              'charge'     : atom.charge   
                              })
            
            representations = {}
            
            for rep_name, rep_data in vobj.representations.items():
                if rep_data == None:
                    pass
                
                else:               
                    if rep_name == 'spheres':
                        representations[rep_name] = list(rep_data.atomic_indexes)
                    else:
                        
                        # this try / except  should be replaced later
                        
                        try:#else:
                            print(vobj.name, vobj.index, rep_name, rep_data.indexes )
                            representations[rep_name] = list(rep_data.indexes)
                        except:
                            
                            #print()
                            pass
            

            #'''
            vismol_objects_dic[vobj_id] = {'name'                 : vobj.name                 ,
                                           'atoms'                : atoms                     ,
                                           'index_bonds'          : list(vobj.index_bonds)    ,
                                           'dynamic_bons'         : list(vobj.dynamic_bons)   ,
                                           'active'               : vobj.active               ,
                                           'frames'               : vobj.frames               ,
                                           'vobj_id'              : vobj_id                   ,
                                           'easyhybrid_system_id' : vobj.easyhybrid_system_id ,
                                           'representations'      : representations           ,
                                           'color_palette'        : vobj.color_palette        ,
                                           }
            
            
            
            
        #********************************************************************** 
        
        
        
        #---------------------------------------------------------------
        #                 P D Y N A M O    P R O J E C T S
        #---------------------------------------------------------------
        pdynamo_projects = {
                            'name'                   : self.pDynamo_session.name,
                            'nbModel_default'        : self.pDynamo_session.nbModel_default,
                            'fixed_color'            : self.pDynamo_session.fixed_color,
                            'pdynamo_distance_safety': self.pDynamo_session.pdynamo_distance_safety,
                            'active_id'              : self.pDynamo_session.active_id,
                            'counter'                : self.pDynamo_session.counter  ,
                            'color_palette_counter'  : self.pDynamo_session.color_palette_counter,
                            'systems'                : {0:None},
                            }
        
        
        for system_id, system in self.pDynamo_session.systems.items():
            pdynamo_projects['systems'][system_id] = {}
            for key in system.keys():
                
                if key == 'vismol_object':
                    pdynamo_projects['systems'][system_id]['vismol_object'] = system['vismol_object'].index
                else:
                    pdynamo_projects['systems'][system_id][key] = system[key]
            
            '''
            pdynamo_projects['systems'][system_id] = {
                                                    'id'            : system['id'           ],
                                                    'name'          : system['name'         ],
                                                    'system'        : system['system'       ],
                                                    'vismol_object' : system['vismol_object'].index,
                                                    'active'        : system['active'       ],
                                                    'bonds'         : system['bonds'        ],
                                                    'sequence'      : system['sequence'     ],
                                                    'qc_table'      : system['qc_table'     ],
                                                    'fixed_table'   : system['fixed_table'  ],
                                                    'color_palette' : system['color_palette'],
                                                    }
            '''
        #---------------------------------------------------------------
        #               V I S M O L   S E S S I O N
        #---------------------------------------------------------------
        
        vm_session = {'vobj_counter': self.vobj_counter}

        easyhybrid_session_data = { 
                                    'pdynamo_projects'   : pdynamo_projects  ,
                                    'vismol_objects_dic' : vismol_objects_dic,
                                    'vm_session'         : vm_session    ,
                                    }
        
        with open(filename,'wb') as outfile:
            pickle.dump(easyhybrid_session_data, outfile)
        



    def load_session (self, infile):
        """ Function doc """
        self.pDynamo_session = self.main_session.pDynamo_session
        infile             = open(infile,'rb')
        easyhybrid_session_data = pickle.load(infile)
        
        '''
        self.main_session.vm_session.vismol_objects = []
        #---------------------------------------------------------------
        print(self.main_session.vm_session.vm_session_vbos)
        for index in self.main_session.vm_session.vm_session_vbos:
            #GL.glDeleteBuffers(1, index)
            idn_array = []
            idn_array.append(index)
            array = (ctypes.c_int * len(idn_array))(*idn_array)
            GL.glDeleteBuffers( 1 , ctypes.byref( array) )
        
        
        #for vismol_object in self.main_session.vm_session.vismol_objects:
            
        #---------------------------------------------------------------
        #'''
        

        
        self.vobj_counter = easyhybrid_session_data['vm_session']['vobj_counter']
        
        self.pDynamo_session.name            = easyhybrid_session_data['pdynamo_projects']['name']
        self.pDynamo_session.nbModel_default = easyhybrid_session_data['pdynamo_projects']['nbModel_default']
        self.pDynamo_session.fixed_color     = easyhybrid_session_data['pdynamo_projects']['fixed_color']
        
        
        self.pDynamo_session.pdynamo_distance_safety = easyhybrid_session_data['pdynamo_projects']['pdynamo_distance_safety']
        self.pDynamo_session.active_id               = easyhybrid_session_data['pdynamo_projects']['active_id']
        self.pDynamo_session.counter                 = easyhybrid_session_data['pdynamo_projects']['counter']
        self.pDynamo_session.systems                 = easyhybrid_session_data['pdynamo_projects']['systems']
        self.pDynamo_session.color_palette_counter   = easyhybrid_session_data['pdynamo_projects']['color_palette_counter']
        
        
        
        
        for vobj_id, vobject_data in easyhybrid_session_data['vismol_objects_dic'].items():

            frames          = vobject_data['frames']
            name            = vobject_data['name']
            bonds           = vobject_data['index_bonds']
            dynamic_bons    = vobject_data['dynamic_bons']
            atoms           = vobject_data['atoms']
            representations = vobject_data['representations']
            color_palette   = vobject_data['color_palette']
            #print('atoms:', len(atoms), 'coords', len(list(frames[0]))/3)
        
            vismol_object  = VismolObject.VismolObject(
                                                       #active                         = easyhybrid_session_data['vismol_objects_dic'][vobj_id]['active']
                                                       name                           = name, 
                                                       atoms                          = atoms, 
                                                       vm_session                  = self, 
                                                       bonds_pair_of_indexes          = list(bonds),
                                                       auto_find_bonded_and_nonbonded = False,
                                                       trajectory                     = frames,
                                                       color_palette                  = color_palette,
                                                       )
                
            vismol_object.index                = vobj_id
            vismol_object.active               = vobject_data['active']
            vismol_object.dynamic_bons         = vobject_data['dynamic_bons']
            vismol_object.easyhybrid_system_id = vobject_data['easyhybrid_system_id']
            vismol_object.set_model_matrix(self.glwidget.vm_widget.model_mat)
            
            # - - - - - - - - - R E P R E S E N T A T I O N - - - - - - - - - - - - - - - 
            #for rep_key in representations.keys():
                
            self.add_vismol_object_to_vismol_session (pdynamo_session = self.pDynamo_session, 
                                                      rep             = representations, 
                                                      vismol_object   = vismol_object, 
                                                      vobj_count      = False,
                                                      autocenter      = True)
        
        
        
        for key, system in self.pDynamo_session.systems.items():
            system['vismol_object'] =   self.vismol_objects_dic[system['vismol_object']]
        

        self.pDynamo_session.refresh_qc_and_fixed_representations() 

        

