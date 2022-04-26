import numpy as np
import math
import ctypes
from OpenGL import GL
import time
import glCore.sphere_data as sphd
import glCore.cylinder_data as cyd
import glCore.matrix_operations as mop
import vModel.cartoon as cartoon
import ctypes

from  vModel.MolecularProperties import residues_dictionary

#from   glCore.sphere_representation import _create_frame_sphere_data

#import EDTSurf.edtsurf as  edtsurf
'''
def _create_frame_sphere_data (frame, atoms ,offset, elems, scale, level,qtty ):
    """ Function doc """
    coords  = sphd.sphere_vertices[level]*qtty
    centers = sphd.sphere_vertices[level]*qtty
    for a, atom in enumerate(atoms):
        pos = atom.coords (frame)
        centers[a*offset:(a+1)*offset] = [pos[0],pos[1],pos[2]]*elems
        
        for i in range(elems):
            coords[a*offset+i*3]   *= atom.radius * scale
            coords[a*offset+i*3+1] *= atom.radius * scale
            coords[a*offset+i*3+2] *= atom.radius * scale
            coords[a*offset+i*3]   += pos[0]
            coords[a*offset+i*3+1] += pos[1]
            coords[a*offset+i*3+2] += pos[2]
    
    coords  = np.array(coords, dtype=np.float32)
    centers = np.array(centers, dtype=np.float32)
    
    return coords, centers
    
    
    #self.centers_list.append(self.centers)
    #self.frames.append(self.coords)
'''










class Representation:
    """ Class doc """

    def _make_gl_VAO (self):
        """ Function doc """
        vao = GL.glGenVertexArrays(1)
        GL.glBindVertexArray(vao)
        return vao

    def _make_gl_index_buffer(self, indexes):
        """ Function doc """
        ind_vbo = GL.glGenBuffers(1)
        self.vbos_list.append(ind_vbo)

        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indexes.nbytes, indexes, GL.GL_DYNAMIC_DRAW)
        return ind_vbo

    def _make_gl_coord_buffer(self, coords, program):
        """ Function doc """
        coord_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.nbytes, coords, GL.GL_STATIC_DRAW)
        att_position = GL.glGetAttribLocation(program, 'vert_coord')
        GL.glEnableVertexAttribArray(att_position)
        GL.glVertexAttribPointer(att_position, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
        self.vbos_list.append(coord_vbo)
        #GL.glDisableVertexAttribArray(att_position)
        return coord_vbo
        
    def _make_gl_normal_buffer(self, normals, program):
        """ Function doc """
        normal_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, normal_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, normals.nbytes, normals, GL.GL_STATIC_DRAW)

        att_normals = GL.glGetAttribLocation(program, 'vert_normal')
        if att_normals > 0:
            GL.glEnableVertexAttribArray(att_normals)
            GL.glVertexAttribPointer(att_normals, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*normals.itemsize, ctypes.c_void_p(0))
        #GL.glDisableVertexAttribArray(att_normals)
        self.vbos_list.append(normal_vbo)
        return normal_vbo        

    def _make_gl_color_buffer(self, colors, program):
        """ Function doc """
        col_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_STATIC_DRAW)
        att_colors = GL.glGetAttribLocation(program, 'vert_color')
        GL.glEnableVertexAttribArray(att_colors)
        GL.glVertexAttribPointer(att_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))
        #GL.glDisableVertexAttribArray(att_colors)
        self.vbos_list.append(col_vbo)
        return col_vbo

    def _make_gl_size_buffer (self, dot_sizes, program):
        """ Function doc """
        size_vbo = GL.glGenBuffers(1)
        self.visObj.vm_session.vm_session_vbos.append(size_vbo)
        
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, size_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, dot_sizes.nbytes, dot_sizes, GL.GL_STATIC_DRAW)
        att_size = GL.glGetAttribLocation(program, 'vert_dot_size')
        GL.glEnableVertexAttribArray(att_size)
        GL.glVertexAttribPointer(att_size, 1, GL.GL_FLOAT, GL.GL_FALSE, dot_sizes.itemsize, ctypes.c_void_p(0))
        #GL.glDisableVertexAttribArray(att_size)
        self.vbos_list.append(size_vbo)
        return size_vbo
    
    
    def _set_colors_to_buffer (self, col_vbo = True):
        """ Function doc """
        #try:
        frame = self.visObj.colors
        #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, visObj.line_buffers[1])
        frame = np.array(self.visObj.colors, dtype = np.float32)
        
        if col_vbo:
                if self.col_vbo:
                    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 
                                    self.col_vbo    )
                    
                    GL.glBufferData(GL.GL_ARRAY_BUFFER, 
                                    frame.nbytes      ,
                                    frame             , 
                                    GL.GL_STATIC_DRAW)   
                    #except:
                    #    
                    #   #print ('wrong type:', self.col_vbo, type(self.col_vbo))
                else: 
                    pass
        else: 
            pass
        #except:
        #   print('_set_colors_to_buffer -  error')
        
        
        
    def _set_coordinates_to_buffer (self, coord_vbo = True, sel_coord_vbo = True):
        ''' This function assigns the coordinates to 
        be drawn by the function  "draw_representation"'''
        
        frame = self.glCore._safe_frame_exchange(self.visObj)
        #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, visObj.line_buffers[1])
        
        if coord_vbo:
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 
                            self.coord_vbo    )
            
            GL.glBufferData(GL.GL_ARRAY_BUFFER, 
                            frame.nbytes      ,
                            frame             , 
                            GL.GL_STATIC_DRAW)   
        else:
            pass
        
        
        
        if sel_coord_vbo:
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 
                            self.sel_coord_vbo )
            
            
            GL.glBufferData(GL.GL_ARRAY_BUFFER, 
                            frame.nbytes      ,
                            frame             , 
                            GL.GL_STATIC_DRAW)    

        else:
            pass
            
    def _make_gl_representation_vao_and_vbos (self, 
                                              indexes    = None,
                                              coords     = None,
                                              colors     = None,
                                              dot_sizes  = None,
                                              normals    = None
                                              ):
        """ Function doc """
       #print ('building', self.name,' VAO  and VBOs')    
        self.vao        =   self._make_gl_VAO()
        self.ind_vbo    =   self._make_gl_index_buffer( indexes                        )
        self.coord_vbo  =   self._make_gl_coord_buffer( coords   , self.shader_program )
        self.col_vbo    =   self._make_gl_color_buffer( colors   , self.shader_program )
        
        if dot_sizes is not None:
            self.dot_sizes    =   self._make_gl_size_buffer(np.array(dot_sizes, dtype = np.float32), self.shader_program )
        
        if normals is not None and self.name == "surface":
            self.norm_vbo   =   self._make_gl_normal_buffer ( normals , self.shader_program )
        else:
            pass
        
    
    def _make_gl_sel_representation_vao_and_vbos (self, 
                                                  indexes    = None,
                                                  coords     = None,
                                                  colors     = None,
                                                  dot_sizes  = None,
                                                  ):
        """ Function doc """
       #print ('building', self.name,'background selection  VAO  and VBOs')    
        self.sel_vao        =   self._make_gl_VAO()
        self.sel_ind_vbo    =   self._make_gl_index_buffer( indexes                             )
        self.sel_coord_vbo  =   self._make_gl_coord_buffer( coords    , self.sel_shader_program )
        self.sel_col_vbo    =   self._make_gl_color_buffer( colors    , self.sel_shader_program )
        if dot_sizes is not None:
            self.sel_size_vbo   =   self._make_gl_size_buffer ( np.array(dot_sizes, dtype = np.float32) , self.sel_shader_program )
        else:
            pass
        
        

    def _check_VAO_and_VBOs (self, indexes = None):
        """ Function doc """
        if self.sel_vao is None:
            self._make_gl_vao_and_vbos ()
        else:
            pass
    
    def _enable_anti_alis_to_lines (self):
        """ Function doc """
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        GL.glEnable(GL.GL_LINE_SMOOTH)
        GL.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST)

    
    def define_new_indexes_to_VBO ( self, input_indexes = []):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        
        
        if input_indexes == []:
            self.active = False
            return None
        else:
            self.active = True
            
        self.indexes = input_indexes
        self.indexes = np.array(self.indexes,dtype=np.uint32)
        
        #if self.sel_ind_vbo:
        #    idn_array = []
        #    idn_array.append(self.sel_ind_vbo   )
        #    array = (ctypes.c_int * len(idn_array))(*idn_array)
        #    GL.glDeleteBuffers( 1 , ctypes.byref( array) )
        #    
        #    self.ind_vbo = GL.glGenBuffers(1)
        #    self.visObj.vm_session.vm_session_vbos.append(self.ind_vbo)
        #
        #    #ind_vbo = self.sel_ind_vbo
        #    GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.sel_ind_vbo)
        #    GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indexes.nbytes, self.indexes, GL.GL_DYNAMIC_DRAW)
        #
        #if self.ind_vbo:
        #    idn_array = []
        #    idn_array.append(self.ind_vbo   )
        #    array = (ctypes.c_int * len(idn_array))(*idn_array)
        #    GL.glDeleteBuffers( 1 , ctypes.byref( array) )
        #
        #    self.ind_vbo = GL.glGenBuffers(1)
        #    self.visObj.vm_session.vm_session_vbos.append(self.ind_vbo)
        #    
        #    GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ind_vbo)
        #    GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indexes.nbytes, self.indexes, GL.GL_DYNAMIC_DRAW)
        
        

        #ind_vbo = self.ind_vbo
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indexes.nbytes, self.indexes, GL.GL_DYNAMIC_DRAW)
        
        #ind_vbo = self.sel_ind_vbo
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.sel_ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indexes.nbytes, self.indexes, GL.GL_DYNAMIC_DRAW)

    def change_vbo_colors  (self,  colors = []):
        """ Function doc """
        colors = np.array(colors,dtype=np.float32)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_STATIC_DRAW)
        att_colors = GL.glGetAttribLocation(self.shader_program, 'vert_color')
        GL.glEnableVertexAttribArray(att_colors)
        GL.glVertexAttribPointer(att_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))


    def delete_buffers (self, ind_vbo = True, coord_vbo = True, vao = True, col_vbo = True, indexes = []):
        for vbo in self.vbos_list:
            try:
                idn_array = []
                idn_array.append(vbo)
                array = (ctypes.c_int * len(idn_array))(*idn_array)
                GL.glDeleteBuffers( 1 , ctypes.byref( array) )
            except:
               print('fail in delete index: ', vbo )

class LinesRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'lines', active = True, _type = 'mol', visObj = None, glCore = None, indexes = []):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type

        self.visObj             = visObj
        self.vm_session      = visObj.vm_session
        
        self.glCore             = glCore
        self.indexes            = indexes
        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        self.vbos_list =[]


    def _make_gl_vao_and_vbos (self, indexes = None):
        """ Function doc """
        
        #if indexes is not None:
        #   #print ('_make_gl_vao_and_vbos',indexes)
        #    pass
        #else:
        #   #print ('_make_gl_vao_and_vbos',indexes)
        #    indexes = np.array(self.visObj.index_bonds,dtype=np.uint32)

        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        
        if indexes:
            self.indexes = np.array(indexes,dtype=np.uint32)
        else:
            self.indexes = np.array(self.visObj.index_bonds,dtype=np.uint32)
        
        #print(self.visObj.index_bonds,self.visObj.frames[0], self.visObj.colors )
        #indexes = np.array(self.visObj.index_bonds,dtype=np.uint32)
        coords  = self.visObj.frames[0]
        colors  = self.visObj.colors
        #colors  = self.visObj.colors_rainbow

        self._make_gl_representation_vao_and_vbos (indexes    = self.indexes,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   )
        colors_idx = self.visObj.color_indexes
        self._make_gl_sel_representation_vao_and_vbos (indexes    = self.indexes    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )



    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        self._enable_anti_alis_to_lines()
        GL.glUseProgram(self.shader_program)
        
        line_width = self.visObj.vm_session.vConfig.gl_parameters['line_width']
        #print('drawing lines')
        line_width = (line_width*200/abs(self.glCore.dist_cam_zrp)/2)**0.5  #40/abs(self.glCore.dist_cam_zrp)
        GL.glLineWidth(line_width)


        self.glCore.load_matrices(self.shader_program, self.visObj.model_mat)
        self.glCore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
            #print(self.visObj.name)
            #self.define_new_indexes_to_VBO ( self.visObj.index_bonds)
            GL.glDrawElements(GL.GL_LINES, int(len(self.visObj.index_bonds)*1), GL.GL_UNSIGNED_INT, None)

        GL.glBindVertexArray(0)
        #GL.glLineWidth(1)
        GL.glUseProgram(0)
        GL.glDisable(GL.GL_LINE_SMOOTH)
        GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self, line_width_factor = 5):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        #if self.sel_vao is None:
        #   #print ('_make_gl_vao_and_vbos')    
        #    self._make_gl_vao_and_vbos ()
        #else:
        #    pass
        line_width = self.visObj.vm_session.vConfig.gl_parameters['line_width_selection'] 
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(line_width)#*line_width_factor) # line_width_factor -> turn the lines thicker

        self.glCore.load_matrices(self.sel_shader_program, self.visObj.model_mat)
        GL.glBindVertexArray(self.sel_vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea
            '''
            self._set_coordinates_to_buffer (coord_vbo = False, sel_coord_vbo = True)

            #frame = self.glCore._safe_frame_exchange(self.visObj)
            #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.sel_coord_vbo)
            #
            #GL.glBufferData(GL.GL_ARRAY_BUFFER, frame.nbytes,
            #                frame, 
            #                GL.GL_STATIC_DRAW)              
            GL.glDrawElements(GL.GL_LINES, int(len(self.visObj.index_bonds)*1), GL.GL_UNSIGNED_INT, None)  
        GL.glBindVertexArray(0)
        GL.glLineWidth(2)
        GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)



        
class DynamicBonds (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'sticks', active = True, _type = 'mol', visObj = None, glCore = None,  indexes = [] ):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type

        self.visObj             = visObj
        self.glCore             = glCore
        self.vbos_list =[]

        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        
        #visObj.find_dynamic_bonds (atom_list = selection.selected_atoms, index_list = None, update = True )
        
        if indexes == []:
            self.indexes = np.array(self.visObj.index_bonds, dtype=np.uint32)
        else:
            self.indexes = np.array(indexes, dtype=np.uint32)
            
    def _make_gl_vao_and_vbos (self, indexes = [], all_white = True):
        """ Function doc """
        #if indexes == []:
        #    self.indexes = np.array(self.visObj.index_bonds, dtype=np.uint32)
        #else:
        #    self.indexes = np.array(indexes, dtype=np.uint32)
        
        #self.indexes = np.array([0,1,0,2,1,2], dtype=np.uint32)
        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        
        
        #indexes = np.array(self.visObj.index_bonds, dtype=np.uint32)
        coords  = self.visObj.frames[0]

        if all_white:
            colors  = np.array( [1.0]*len(self.visObj.colors), dtype=np.float32)
        else:
            colors  = self.visObj.colors

        self._make_gl_representation_vao_and_vbos (indexes    = self.indexes,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   )
        colors_idx = self.visObj.color_indexes
        self._make_gl_sel_representation_vao_and_vbos (indexes    = self.indexes    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        self._enable_anti_alis_to_lines()

        GL.glUseProgram(self.shader_program)
        GL.glLineWidth(40/abs(self.glCore.dist_cam_zrp))

        self.glCore.load_matrices(self.shader_program, self.visObj.model_mat)
        self.glCore.load_fog(self.shader_program)
        self.glCore.load_lights(self.shader_program)

        GL.glBindVertexArray(self.vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            
            frame = self.glCore.frame
            #try:
            #print (frame, self.visObj.dynamic_bonds[frame])
            #self.define_new_indexes_to_VBO ( self.visObj.index_bonds)
            if frame < len(self.visObj.dynamic_bonds):
                self.define_new_indexes_to_VBO ( self.visObj.dynamic_bonds[frame])
                self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
                #print("line 515", self.visObj.dynamic_bonds[frame])
                GL.glDrawElements(GL.GL_LINES, int(len(self.visObj.dynamic_bonds[frame])*1), GL.GL_UNSIGNED_INT, None)
            else:
                self.define_new_indexes_to_VBO ( self.visObj.dynamic_bonds[-1])
                self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
                #print("line 520")
                GL.glDrawElements(GL.GL_LINES, int(len(self.visObj.dynamic_bonds[-1])*1), GL.GL_UNSIGNED_INT, None)
            #except:
            #    pass
            

        
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()

        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(20)
        GL.glDisable(GL.GL_LINE_SMOOTH)
        GL.glDisable(GL.GL_BLEND)

        self.glCore.load_matrices(self.sel_shader_program, self.visObj.model_mat)
        GL.glBindVertexArray(self.sel_vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea
            '''
            frame = self.glCore.frame

            if frame < len(self.visObj.dynamic_bonds):
                self.define_new_indexes_to_VBO ( self.visObj.dynamic_bonds[frame])
                self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
                GL.glDrawElements(GL.GL_LINES, int(len(self.visObj.dynamic_bonds[frame])), GL.GL_UNSIGNED_INT, None)
            else:
                self.define_new_indexes_to_VBO ( self.visObj.dynamic_bonds[-1])
                self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
                GL.glDrawElements(GL.GL_LINES, int(len(self.visObj.dynamic_bonds[-1])), GL.GL_UNSIGNED_INT, None)
            
            
            
            #frame = self.glCore.frame
            
            
            ##try:
            ##print (frame, self.visObj.dynamic_bonds[frame])
            ##self.define_new_indexes_to_VBO ( self.visObj.index_bonds)
            #self.define_new_indexes_to_VBO ( self.visObj.dynamic_bonds[frame])
            #self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
            #GL.glDrawElements(GL.GL_LINES, int(len(self.visObj.dynamic_bonds[frame])*2), GL.GL_UNSIGNED_INT, None)
            
            #self._set_coordinates_to_buffer (coord_vbo = False, sel_coord_vbo = True)
            #GL.glDrawElements(GL.GL_LINES, int(len(self.visObj.index_bonds)*2), GL.GL_UNSIGNED_INT, None)
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)
        GL.glDisable(GL.GL_DEPTH_TEST)


        
        
        
class SticksRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'sticks', active = True, _type = 'mol', visObj = None, glCore = None,  indexes = [] ):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type

        self.visObj             = visObj
        self.glCore             = glCore
        self.vbos_list =[]

        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        
        if indexes == [] or indexes == None:
            self.indexes = np.array([], dtype=np.uint32)
        else:
            self.indexes = np.array(indexes, dtype=np.uint32)

    def delete_buffers (self, ind_vbo = True, coord_vbo = True, vao = True, col_vbo = True, indexes = []):
        """ Function doc """
        
        delete_buffers = {
                           'vao'          : self.vao           ,
                           'ind_vbo'      : self.ind_vbo       ,
                           'coord_vbo'    : self.coord_vbo     ,
                           'col_vbo'      : self.col_vbo       ,
                           'size_vbo'     : self.size_vbo      ,
                           'sel_vao'      : self.sel_vao       ,
                           'sel_ind_vbo'  : self.sel_ind_vbo   ,
                           'sel_coord_vbo': self.sel_coord_vbo ,
                           'sel_col_vbo'  : self.sel_col_vbo   ,
                           'sel_size_vbo' : self.sel_size_vbo  }
        
        #
        #for key , index in delete_buffers.items():
        #    try:
        #        idn_array = []
        #        idn_array.append(index)
        #        array = (ctypes.c_int * len(idn_array))(*idn_array)
        #        GL.glDeleteBuffers( 1 , ctypes.byref( array) )
        #    except:
        #       #print('fail in delete index: ', index, key)
                
        #idn_array = []
        #idn_array.append(self.vao   )
        #array = (ctypes.c_int * len(idn_array))(*idn_array)
        #GL.glDeleteBuffers( 1 , ctypes.byref( array) )

        
        
        
        #idn_array = []
        #idn_array.append(self.sel_coord_vbo   )
        #array = (ctypes.c_int * len(idn_array))(*idn_array)
        #GL.glDeleteBuffers( 1 , ctypes.byref( array) )
        #
        #idn_array = []
        #idn_array.append(self.sel_col_vbo   )
        #array = (ctypes.c_int * len(idn_array))(*idn_array)
        #GL.glDeleteBuffers( 1 , ctypes.byref( array) )
        #
        #
        #idn_array = []
        #idn_array.append(self.vao   )
        #array = (ctypes.c_int * len(idn_array))(*idn_array)
        #GL.glDeleteBuffers( 1 , ctypes.byref( array) )
        '''
        idn_array = []
        idn_array.append(self.ind_vbo   )
        array = (ctypes.c_int * len(idn_array))(*idn_array)
        GL.glDeleteBuffers( 1 , ctypes.byref( array) )
        '''
        
        '''
        idn_array = []
        idn_array.append(self.coord_vbo   )
        array = (ctypes.c_int * len(idn_array))(*idn_array)
        GL.glDeleteBuffers( 1 , ctypes.byref( array) )
        '''
        
        '''
        idn_array = []
        idn_array.append(self.col_vbo   )
        array = (ctypes.c_int * len(idn_array))(*idn_array)
        GL.glDeleteBuffers( 1 , ctypes.byref( array) )
        '''
    
    def _make_gl_vao_and_vbos (self, indexes = []):
        """ Function doc """
        #if indexes == []:
        #    self.indexes = np.array(self.visObj.index_bonds, dtype=np.uint32)
        #else:
        #    self.indexes = np.array(indexes, dtype=np.uint32)
        
        #self.indexes = np.array([0,1,0,2,1,2], dtype=np.uint32)
        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        

        #indexes = np.array(self.visObj.index_bonds, dtype=np.uint32)
        coords  = self.visObj.frames[0]
        colors  = self.visObj.colors

        self._make_gl_representation_vao_and_vbos (indexes    = self.indexes,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   )
        colors_idx = self.visObj.color_indexes
        self._make_gl_sel_representation_vao_and_vbos (indexes    = self.indexes    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        self._enable_anti_alis_to_lines()

        GL.glUseProgram(self.shader_program)
        GL.glLineWidth(40/abs(self.glCore.dist_cam_zrp))
        #print('drawing sticks')

        self.glCore.load_matrices(self.shader_program, self.visObj.model_mat)
        self.glCore.load_fog(self.shader_program)
        self.glCore.load_lights(self.shader_program)

        GL.glBindVertexArray(self.vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
            try:
                GL.glDrawElements(GL.GL_LINES, len(self.visObj.index_bonds), GL.GL_UNSIGNED_INT, None)
            except:
                pass
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()

        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(100)
        GL.glDisable(GL.GL_LINE_SMOOTH)
        GL.glDisable(GL.GL_BLEND)

        self.glCore.load_matrices(self.sel_shader_program, self.visObj.model_mat)
        GL.glBindVertexArray(self.sel_vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea
            '''
            self._set_coordinates_to_buffer (coord_vbo = False, sel_coord_vbo = True)
            try:
                GL.glDrawElements(GL.GL_LINES, len(self.visObj.index_bonds), GL.GL_UNSIGNED_INT, None)
            except:
                pass
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)
        GL.glDisable(GL.GL_DEPTH_TEST)




class RibbonsRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'ribbon', active = True, _type = 'mol', visObj = None, glCore = None):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type

        self.visObj             = visObj
        self.glCore             = glCore
        self.vbos_list =[]

        
        if self.visObj.c_alpha_bonds == []:
            self.visObj.get_backbone_indexes ()
            
            if self.visObj.c_alpha_bonds == []:
                self.active =  False
        
        else:
            pass
        
        
        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        #print('opaaaa')

    def _make_gl_vao_and_vbos (self, indexes = None):
        """ Function doc """
        
        #if indexes is not None:
        #   #print ('_make_gl_vao_and_vbos',indexes)
        #    pass
        #else:
        #   #print ('_make_gl_vao_and_vbos',indexes)
        #    indexes = np.array(self.visObj.index_bonds,dtype=np.uint32)

        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        
        indexes = []
        
        # c_alpha_bonds
        for bond in self.visObj.c_alpha_bonds:
            i = bond.atom_index_i
            j = bond.atom_index_j
            #print (i, j)
            indexes.append(i)
            indexes.append(j)
            
        if indexes == []:
            self.activate = False
        else:
            
            indexes = np.array(indexes,dtype=np.uint32)

            coords  = self.visObj.frames[0]
            #colors  = self.visObj.colors
            colors  = self.visObj.colors_rainbow
            #colors  = np.array([1.0 ]*len(coords),dtype=np.float32)
            self._make_gl_representation_vao_and_vbos (indexes    = indexes,
                                                       coords     = coords ,
                                                       colors     = colors ,
                                                       dot_sizes  = None   ,
                                                       )
            colors_idx = self.visObj.color_indexes
            #colors_idx = np.array([1.0, 1.0 , 1.0 ],dtype=np.float32)
            self._make_gl_sel_representation_vao_and_vbos (indexes    = indexes    ,
                                                           coords     = coords     ,
                                                           colors     = colors_idx ,
                                                           dot_sizes  = None       ,
                                                           )



    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        self._enable_anti_alis_to_lines()
        GL.glUseProgram(self.shader_program)
        
        ribbon_width = self.visObj.vm_session.vConfig.gl_parameters['ribbon_width']
        LineWidth = ((ribbon_width*80)/abs(self.glCore.dist_cam_zrp)/2)  #40/abs(self.glCore.dist_cam_zrp)
        #print(LineWidth)
        GL.glLineWidth(LineWidth)


        self.glCore.load_matrices(self.shader_program, self.visObj.model_mat)
        self.glCore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
            GL.glDrawElements(GL.GL_LINES, int(len(self.visObj.index_bonds)), GL.GL_UNSIGNED_INT, None)

        GL.glBindVertexArray(0)
        #GL.glLineWidth(1)
        GL.glUseProgram(0)
        GL.glDisable(GL.GL_LINE_SMOOTH)
        GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        #if self.sel_vao is None:
        #   #print ('_make_gl_vao_and_vbos')    
        #    self._make_gl_vao_and_vbos ()
        #else:
        #    pass
        line_width = self.visObj.vm_session.vConfig.gl_parameters['line_width_selection'] 
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(line_width)

        self.glCore.load_matrices(self.sel_shader_program, self.visObj.model_mat)
        GL.glBindVertexArray(self.sel_vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea
            '''
            self._set_coordinates_to_buffer (coord_vbo = False, sel_coord_vbo = True)

            #frame = self.glCore._safe_frame_exchange(self.visObj)
            #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.sel_coord_vbo)
            #
            #GL.glBufferData(GL.GL_ARRAY_BUFFER, frame.nbytes,
            #                frame, 
            #                GL.GL_STATIC_DRAW)              
            #print('aquioh')
            GL.glDrawElements(GL.GL_LINES, int(len(self.visObj.index_bonds)), GL.GL_UNSIGNED_INT, None)  
        GL.glBindVertexArray(0)
        GL.glLineWidth(1)
        GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)




class NonBondedRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'nonbonded', active = True, _type = 'mol', visObj = None, glCore = None, indexes = None):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type
        self.vbos_list =[]

        self.visObj             = visObj
        self.glCore             = glCore
        self.indexes            = indexes
              
        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None


    def _make_gl_vao_and_vbos (self, indexes = None):
        """ Function doc """
        if indexes is not None:
            self.indexes =  indexes
            pass        
        
        else:
            self.indexes = np.array(self.visObj.non_bonded_atoms, dtype=np.uint32)
        
        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        
        #indexes = np.array(self.visObj.index_bonds, dtype=np.uint32)
        #indexes = np.array(self.visObj.non_bonded_atoms, dtype=np.uint32)
        coords  = self.visObj.frames[0]
        colors  = self.visObj.colors

        self._make_gl_representation_vao_and_vbos (indexes    = self.indexes,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   )
        colors_idx = self.visObj.color_indexes
        self._make_gl_sel_representation_vao_and_vbos (indexes    = self.indexes    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        self._enable_anti_alis_to_lines()

        line_width = self.visObj.vm_session.vConfig.gl_parameters['line_width']
        
        GL.glUseProgram(self.shader_program)
        GL.glLineWidth(line_width*20/abs(self.glCore.dist_cam_zrp))

        self.glCore.load_matrices(self.shader_program, self.visObj.model_mat)
        self.glCore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
            GL.glDrawElements(GL.GL_POINTS, int(len(self.visObj.non_bonded_atoms)), GL.GL_UNSIGNED_INT, None)

        GL.glBindVertexArray(0)
        GL.glLineWidth(1)
        GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self, line_width_factor = 5):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(20)

        self.glCore.load_matrices(self.sel_shader_program, self.visObj.model_mat)
        GL.glBindVertexArray(self.sel_vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea
            '''
            self._set_coordinates_to_buffer (coord_vbo = False, sel_coord_vbo = True)
            GL.glDrawElements(GL.GL_POINTS, int(len(self.visObj.non_bonded_atoms)), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glLineWidth(1)
        GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)



class DotsRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'dots', active = True, _type = 'mol', visObj = None, glCore = None, indexes = []):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type

        self.visObj             = visObj
        self.glCore             = glCore
        
        self.vbos_list =[]


        
        
        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        
        if indexes == []:
            self.indexes = np.array(range(len(self.visObj.atoms)), dtype=np.uint32)
        else:
            self.indexes = np.array(indexes, dtype=np.uint32)

    def _make_gl_vao_and_vbos (self, indexes = None):
        """ Function doc """
        #if indexes is not None:
        #    pass
        #else:
        
        #dot_qtty  = int(len(self.visObj.frames[0])/3)
        #indexes = []
        #for i in range(dot_qtty):
        #    indexes.append(i)
        
        indexes = self.indexes       
        
        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        
        #indexes = np.array(self.visObj.index_bonds, dtype=np.uint32)
        #indexes = np.array(self.visObj.idex, dtype=np.uint32)
        coords  = self.visObj.frames[0]
        colors  = self.visObj.colors
        radiues = self.visObj.cov_radiues_list

        self._make_gl_representation_vao_and_vbos (indexes    = indexes,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = radiues,
                                                   )
        colors_idx = self.visObj.color_indexes
        self._make_gl_sel_representation_vao_and_vbos (indexes    = indexes    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = radiues    ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        self._enable_anti_alis_to_lines()
        #print ('DotsRepresentation')
        height = self.visObj.vm_session.glwidget.vm_widget.height

        GL.glUseProgram(self.shader_program)
        #1*self.height dot_size
        #GL.glLineWidth(40/abs(self.glCore.dist_cam_zrp))
        #GL.glPointSize(1.5*height/abs(self.glCore.dist_cam_zrp)) # dot size not included yet
        self.glCore.load_matrices(self.shader_program, self.visObj.model_mat)
        self.glCore.load_lights(self.shader_program)
        
        xyz_coords = self.glCore.glcamera.get_modelview_position(self.visObj.model_mat)
        u_campos   = GL.glGetUniformLocation(self.shader_program, "u_campos")
        GL.glUniform3fv(u_campos, 1, xyz_coords)
        
        
        self.glCore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            #print('aqui')
            self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
            GL.glDrawElements(GL.GL_POINTS, int(len(self.visObj.atoms)), GL.GL_UNSIGNED_INT, None)

        #GL.glBindVertexArray(0)
        #GL.glLineWidth(1)
        #GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glUseProgram(self.sel_shader_program)
        GL.glPointSize(200/abs(self.glCore.dist_cam_zrp))
        #GL.glLineWidth(20)
        self.glCore.load_matrices(self.sel_shader_program, self.visObj.model_mat)
        GL.glBindVertexArray(self.sel_vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea
            '''
            #print(self.name,'draw_background_sel_representation')
            self._set_coordinates_to_buffer (coord_vbo = False, sel_coord_vbo = True)
            GL.glDrawElements(GL.GL_POINTS, int(len(self.visObj.atoms)), GL.GL_UNSIGNED_INT, None)



#class SpheresRepresentationOLD (Representation):
#    """ Class doc """
#    
#    def __init__ (self, name = 'spheres', 
#                      active = True, 
#                       _type = 'mol', 
#                      visObj = None, 
#                      glCore = None,  
#                     indexes = None
#                       #atoms = None,
#                       #level = 'level_1', 
#                       ):
#        
#        """ Class initialiser """
#        self.name               = name
#        self.active             = active
#        self.type               = _type
#
#        self.visObj             = visObj
#        self.glCore             = glCore
#        
#        self.atomic_indexes     = indexes
#        
#        
#        # --------------------------------
#        #self.level              = level
#        self.level              = self.visObj.vm_session.vConfig.gl_parameters['sphere_quality']
#        self.scale              = self.visObj.vm_session.vConfig.gl_parameters['sphere_scale']
#        
#        
#        if self.atomic_indexes is None:
#            self.atoms          = self.visObj.atoms
#        else:
#            self.atoms = []
#            for index in self.atomic_indexes:
#                self.atoms.append(self.visObj.atoms[index])
#            
#            #self.atoms          = atoms
#        
#        
#        self.coords             = None
#        self.colors             = None
#        self.centers            = None
#        self.indexes            = None
#        self.triangles          = None
#        self.frames             = []
#        # --------------------------------
#
#        # representation 	
#        self.vao            = None
#        self.ind_vbo        = None
#        self.coord_vbo      = None
#        self.col_vbo        = None
#        self.size_vbo       = None
#           
#
#        # bgrd selection   
#        self.sel_vao        = None
#        self.sel_ind_vbo    = None
#        self.sel_coord_vbo  = None
#        self.sel_col_vbo    = None
#        self.sel_size_vbo   = None
#
#
#        #     S H A D E R S
#        self.shader_program     = None
#        self.sel_shader_program = None
#        #self._check_VAO_and_VBOs()
#    
#    def update_atomic_indexes (self, indexes = None):
#        """ Function doc """
#        self.atomic_indexes = []
#        self.atoms =[]
#        for atom in self.visObj.atoms:
#            if atom.spheres:
#                #print (atom.name ,atom.index,  atom.spheres)
#                index  = atom.index -1
#                self.atomic_indexes.append(index)
#                self.atoms.append(atom)
#        #for index in indexes:
#        #    self.atoms.append(self.visObj.atoms[index])
#        #    self.atomic_indexes.append(index)
#        
#        self._create_sphere_data() 
#        self._update_sphere_data_to_VBOs () 
#        self.active = True
#        #print('1232 update_atomic_indexes')
#        #print(self.atomic_indexes)
#        
#    def _update_sphere_data_to_VBOs (self):
#        """ Function doc """
#        
#        #GL.glDeleteVertexArrays( 1, self.vao)        
#        #GL.glDeleteBuffers(1, self.ind_vbo)
#        #GL.glDeleteBuffers(1, self.coord_vbo)
#        #GL.glDeleteBuffers(1, self.centr_vbo)
#        #GL.glDeleteBuffers(1, self.col_vbo)
#        #self._make_gl_vao_and_vbos ()
#       #print ('self.ind_vbo', self.ind_vbo, self.visObj.name, self.visObj.index )
#        #self._check_VAO_and_VBOs()
#        #if self.ind_vbo is None:
#        #    self.draw_representation()
#        #    self._make_gl_vao_and_vbos()
#        
#        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ind_vbo)
#        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indexes.nbytes, self.indexes, GL.GL_DYNAMIC_DRAW)
#        
#        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.coord_vbo)
#        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.coords.nbytes, self.coords, GL.GL_STATIC_DRAW)
#        
#        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.centr_vbo)
#        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.centers.itemsize*len(self.centers), self.centers, GL.GL_STATIC_DRAW)
#        
#        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
#        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.colors.itemsize*len(self.colors), self.colors, GL.GL_STATIC_DRAW)
#
#    
#    def _create_sphere_data(self ):
#        #print('1258_create_sphere_data')
#        """ Function doc """
#        init = time.time()
#        #cdef Py_ssize_t a, i, qtty, elems, offset, inds_e
#        qtty = int(len( self.atoms))
#        nucleus = [0.0, 0.0, 0.0]*qtty
#        colores = [0.0, 0.0, 0.0]*qtty
#        coords  = sphd.sphere_vertices[self.level]*qtty
#        centers = sphd.sphere_vertices[self.level]*qtty
#        colors  = sphd.sphere_vertices[self.level]*qtty
#        indexes = np.array(sphd.sphere_triangles[self.level]*qtty, dtype=np.uint32)
#        elems  = int(len(sphd.sphere_vertices[self.level])/3)
#        offset = int(len(sphd.sphere_vertices[self.level]))
#        inds_e = int(len(sphd.sphere_triangles[self.level]))
#        
#        #print(coords)
#        self.centers_list = []
#        self.frames = []
#        frame =0
#        #for frame in range(len(self.visObj.frames)-1):
#        
#        for a, atom in enumerate( self.atoms ):
#            pos = atom.coords (frame)
#            #print (pos, atom.index, frame)
#
#            colors[a*offset:(a+1)*offset]  = [atom.color[0],atom.color[1],atom.color[2]]*elems
#            centers[a*offset:(a+1)*offset] = [pos[0],pos[1],pos[2]]*elems
#
#            for i in range(elems):
#                coords[a*offset+i*3]   *= atom.radius * self.scale
#                coords[a*offset+i*3+1] *= atom.radius * self.scale
#                coords[a*offset+i*3+2] *= atom.radius * self.scale
#                coords[a*offset+i*3]   += pos[0]
#                coords[a*offset+i*3+1] += pos[1]
#                coords[a*offset+i*3+2] += pos[2]
#            indexes[a*inds_e:(a+1)*inds_e] += a*elems
#        end = time.time()
#       #print('Time used creating nucleus, vertices and colors:', end-init)
#
#
#
#        self.coords  = np.array(coords, dtype=np.float32)
#        self.frames.append(self.coords)
#        
#        self.centers = np.array(centers, dtype=np.float32)
#        self.centers_list.append(self.centers)
#        self.colors  = np.array(colors, dtype=np.float32)
#        self.indexes = indexes
#        
#        self.triangles = int(len(self.indexes))
#
#        init = time.time()
#
#        if len(self.visObj.frames) > 1:
#            #'''
#            for frame in range(1,len(self.visObj.frames)-1):
#                coords  = sphd.sphere_vertices[self.level]*qtty
#                centers = sphd.sphere_vertices[self.level]*qtty
#                for a, atom in enumerate( self.atoms ):
#                    pos = atom.coords (frame)
#                    centers[a*offset:(a+1)*offset] = [pos[0],pos[1],pos[2]]*elems
#                    
#                    for i in range(elems):
#                        coords[a*offset+i*3]   *= atom.radius * self.scale
#                        coords[a*offset+i*3+1] *= atom.radius * self.scale
#                        coords[a*offset+i*3+2] *= atom.radius * self.scale
#                        coords[a*offset+i*3]   += pos[0]
#                        coords[a*offset+i*3+1] += pos[1]
#                        coords[a*offset+i*3+2] += pos[2]
#                self.coords  = np.array(coords, dtype=np.float32)
#                self.frames.append(self.coords)
#                self.centers = np.array(centers, dtype=np.float32)
#                self.centers_list.append(self.centers)
#            #'''
#        
#        #print (self.centers_list)
#        #print (len(self.centers_list))
#        #print (self.centers)
#        #print (len(self.centers))
#        #print (self.coords)
#        #print (len(self.coords))
#        #print (self.indexes)
#        
#        #coords_menos_centers = []
#        #for i in range(0, len(self.coords)):
#        #    coords_menos_centers.append(self.coords[i] -self.centers[i])
#        #print (coords_menos_centers)
#        
#        
#        end = time.time()
#       #print('Time used creating nucleus, vertices and colors:', end-init)
#    
#        return True
#    
#    '''
#    def _create_sel_sphere_data(self, level):
#        """ Function doc """
#        init = time.time()
#        #cdef Py_ssize_t a, i, qtty, elems, offset, inds_e
#        qtty = int(len(self.visObj.atoms))
#        nucleus = [0.0, 0.0, 0.0]*qtty
#        colores = [0.0, 0.0, 0.0]*qtty
#        coords = sphd.sphere_vertices[level]*qtty
#        colors = sphd.sphere_vertices[level]*qtty
#        indexes = np.array(sphd.sphere_triangles[level]*qtty, dtype=np.uint32)
#        elems = int(len(sphd.sphere_vertices[level])/3)
#        offset = int(len(sphd.sphere_vertices[level]))
#        inds_e = int(len(sphd.sphere_triangles[level]))
#        for a,atom in enumerate(self.visObj.atoms):
#            colors[a*offset:(a+1)*offset] = [atom.color_id[0],atom.color_id[1],atom.color_id[2]]*elems
#            for i in range(elems):
#                coords[a*offset+i*3] *= atom.radius * self.scale
#                coords[a*offset+i*3+1] *= atom.radius * self.scale
#                coords[a*offset+i*3+2] *= atom.radius * self.scale
#                coords[a*offset+i*3] += atom.pos[0]
#                coords[a*offset+i*3+1] += atom.pos[1]
#                coords[a*offset+i*3+2] += atom.pos[2]
#            indexes[a*inds_e:(a+1)*inds_e] += a*elems
#        end = time.time()
#       #print('Time used creating nucleus, vertices and colors for selection:', end-init)
#        self.sel_coords = np.array(coords, dtype=np.float32)
#        self.sel_colors = np.array(colors, dtype=np.float32)
#        self.sel_indexes = indexes
#        return True
#        '''
#    '''
#    def _make_gl_spheres(self, program):
#        """ Function doc """
#        vertex_array_object = GL.glGenVertexArrays(1)
#        GL.glBindVertexArray(vertex_array_object)
#        
#        ind_vbo = GL.glGenBuffers(1)
#        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ind_vbo)
#        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indexes.itemsize*int(len(self.indexes)), self.indexes, GL.GL_DYNAMIC_DRAW)
#        
#        coord_vbo = GL.glGenBuffers(1)
#        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
#        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.coords.itemsize*len(self.coords), self.coords, GL.GL_STATIC_DRAW)
#        gl_coord = GL.glGetAttribLocation(program, 'vert_coord')
#        GL.glEnableVertexAttribArray(gl_coord)
#        GL.glVertexAttribPointer(gl_coord, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*self.coords.itemsize, ctypes.c_void_p(0))
#        
#        centr_vbo = GL.glGenBuffers(1)
#        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, centr_vbo)
#        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.centers.itemsize*len(self.centers), self.centers, GL.GL_STATIC_DRAW)
#        gl_center = GL.glGetAttribLocation(program, 'vert_centr')
#        GL.glEnableVertexAttribArray(gl_center)
#        GL.glVertexAttribPointer(gl_center, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*self.centers.itemsize, ctypes.c_void_p(0))
#        
#        col_vbo = GL.glGenBuffers(1)
#        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, col_vbo)
#        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.colors.itemsize*len(self.colors), self.colors, GL.GL_STATIC_DRAW)
#        gl_colors = GL.glGetAttribLocation(program, 'vert_color')
#        GL.glEnableVertexAttribArray(gl_colors)
#        GL.glVertexAttribPointer(gl_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*self.colors.itemsize, ctypes.c_void_p(0))
#        
#        GL.glBindVertexArray(0)
#        GL.glDisableVertexAttribArray(gl_coord)
#        GL.glDisableVertexAttribArray(gl_center)
#        GL.glDisableVertexAttribArray(gl_colors)
#        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
#        self.spheres_vao = vertex_array_object
#        self.spheres_buffers = (ind_vbo, coord_vbo, col_vbo)
#        self.triangles = int(len(self.indexes))
#        return True
#        '''
#    '''
#    def _make_sel_gl_spheres(self, program):
#        """ Function doc """
#        vertex_array_object = GL.glGenVertexArrays(1)
#        GL.glBindVertexArray(vertex_array_object)
#        
#        ind_vbo = GL.glGenBuffers(1)
#        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ind_vbo)
#        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.sel_indexes.itemsize*int(len(self.sel_indexes)), self.sel_indexes, GL.GL_DYNAMIC_DRAW)
#        
#        coord_vbo = GL.glGenBuffers(1)
#        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
#        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.sel_coords.itemsize*len(self.sel_coords), self.sel_coords, GL.GL_STATIC_DRAW)
#        gl_coord = GL.glGetAttribLocation(program, 'vert_coord')
#        GL.glEnableVertexAttribArray(gl_coord)
#        GL.glVertexAttribPointer(gl_coord, 3, GL.GL_FLOAT, GL.GL_FALSE, self.sel_coords.nbytes, ctypes.c_void_p(0))
#        
#        col_vbo = GL.glGenBuffers(1)
#        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, col_vbo)
#        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.sel_colors.itemsize*len(self.sel_colors), self.sel_colors, GL.GL_STATIC_DRAW)
#        gl_colors = GL.glGetAttribLocation(program, 'vert_color')
#        GL.glEnableVertexAttribArray(gl_colors)
#        GL.glVertexAttribPointer(gl_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, self.sel_colors.nbytes, ctypes.c_void_p(0))
#        
#        GL.glBindVertexArray(0)
#        GL.glDisableVertexAttribArray(gl_coord)
#        GL.glDisableVertexAttribArray(gl_colors)
#        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
#        self.sel_spheres_vao = vertex_array_object
#        self.sel_spheres_buffers = (ind_vbo, coord_vbo, col_vbo)
#        self.sel_triangles = int(len(self.indexes))
#        return True
#    
#
#    '''
#    def _make_gl_vao_and_vbos (self):
#        """ Function doc """
#        #self._create_sphere_data()
#
#        self.shader_program     = self.glCore.shader_programs[self.name]
#        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
#        
#        #indexes = np.array(self.visObj.index_bonds,dtype=np.uint32)
#        #coords  = self.visObj.frames[0]
#        #colors  = self.visObj.colors
#
#        #self._make_gl_representation_vao_and_vbos (indexes    = indexes,
#        #                                           coords     = coords ,
#        #                                           colors     = colors ,
#        #                                           dot_sizes  = None   ,
#        #                                           )
#
#        self.vao = GL.glGenVertexArrays(1)
#        GL.glBindVertexArray(self.vao)
#
#        self.ind_vbo = GL.glGenBuffers(1)
#        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ind_vbo)
#        # RAM -> GPU
#        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indexes.nbytes, self.indexes, GL.GL_DYNAMIC_DRAW)
#        
#        # glDeleteBuffers(self.ind_vbo)
#        #print (self.coords)
#        self.coord_vbo = GL.glGenBuffers(1)
#        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.coord_vbo)
#        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.coords.nbytes, self.coords, GL.GL_STATIC_DRAW)
#        gl_coord = GL.glGetAttribLocation(self.shader_program , 'vert_coord')
#        GL.glEnableVertexAttribArray(gl_coord)
#        GL.glVertexAttribPointer(gl_coord, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*self.coords.itemsize, ctypes.c_void_p(0))
#        
#        self.centr_vbo = GL.glGenBuffers(1)
#        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.centr_vbo)
#        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.centers.itemsize*len(self.centers), self.centers, GL.GL_STATIC_DRAW)
#        gl_center = GL.glGetAttribLocation(self.shader_program , 'vert_centr')
#        GL.glEnableVertexAttribArray(gl_center)
#        GL.glVertexAttribPointer(gl_center, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*self.centers.itemsize, ctypes.c_void_p(0))
#
#        self.col_vbo = GL.glGenBuffers(1)
#        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
#        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.colors.itemsize*len(self.colors), self.colors, GL.GL_STATIC_DRAW)
#        gl_colors = GL.glGetAttribLocation(self.shader_program, 'vert_color')
#        GL.glEnableVertexAttribArray(gl_colors)
#        GL.glVertexAttribPointer(gl_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*self.colors.itemsize, ctypes.c_void_p(0))
#        
#        self.triangles = int(len(self.indexes))
#        
#        
#        
#        colors_idx = self.visObj.color_indexes
#        self._make_gl_sel_representation_vao_and_vbos (indexes    = self.indexes    ,
#                                                       coords     = self.coords     ,
#                                                       colors     = colors_idx ,
#                                                       dot_sizes  = None       ,
#                                                       )
#
#    def draw_representation (self):
#        """ Function doc """
#        self._check_VAO_and_VBOs ()
#        
#        GL.glEnable(GL.GL_DEPTH_TEST)
#        GL.glEnable(GL.GL_CULL_FACE)
#        GL.glCullFace(GL.GL_BACK)
#        
#        GL.glUseProgram          (self.shader_program )
#        self.glCore.load_matrices(self.shader_program , self.visObj.model_mat)
#        self.glCore.load_lights  (self.shader_program )
#        self.glCore.load_fog     (self.shader_program )
#        
#        if self.vao is not None:
#            GL.glBindVertexArray (self.vao)
#            if self.glCore.modified_view:
#                pass
#            
#            else:
#                #self.centers = self.glCore._safe_frame_exchange(self.visObj)
#                '''
#                This function checks if the number of the called frame will not exceed 
#                the limit of frames that each object has. Allowing two objects with 
#                different trajectory sizes to be manipulated at the same time within the 
#                glArea'''
#                #self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
#                frame = self.glCore._get_visObj_frame (self.visObj)
#                if frame >= len(self.frames):
#                    frame = len(self.frames)-1
#                #print(frame)
#                
#                self.coords = self.frames[frame]
#                self.centers =self.centers_list[frame]
#                
#                GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.coord_vbo)
#                GL.glBufferData(GL.GL_ARRAY_BUFFER, self.coords.itemsize*int(len(self.coords)), 
#                                                    self.coords, GL.GL_STATIC_DRAW)
#                
#
#                GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.centr_vbo)
#                GL.glBufferData(GL.GL_ARRAY_BUFFER, self.centers.itemsize*len(self.centers), self.centers, GL.GL_STATIC_DRAW)
#
#
#
#                GL.glDrawElements(GL.GL_TRIANGLES,  self.triangles , GL.GL_UNSIGNED_INT, None)
#        GL.glBindVertexArray(0)
#        GL.glUseProgram(0)
#        GL.glDisable(GL.GL_DEPTH_TEST)
#        
#            
#    def draw_background_sel_representation  (self):
#        """ Function doc """
#        pass
#
#
#
class SpheresRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'Spheres', active = True, _type = 'mol', visObj = None, glCore = None, indexes = None):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type

        self.visObj             = visObj
        self.glCore             = glCore
        self.indexes            = indexes
        self.vbos_list =[]

        self.scale              = 0.07
        
        #self.light_position = np.array([-2.5, 2.5, 3.0],dtype=np.float32)
        #self.light_color = np.array([1.0, 1.0, 1.0, 1.0],dtype=np.float32)
        #self.light_ambient_coef = 0.4
        #self.light_shininess = 5.5
        #self.light_intensity = np.array([0.6, 0.6, 0.6],dtype=np.float32)
        #self.light_specular_color = np.array([1.0, 1.0, 1.0],dtype=np.float32)

        #  'light_position'             : [-2.5, -2.5, 3.0  ] ,
        #  'light_color'                : [ 1.0, 1.0, 1.0,1.0] ,
        #  'light_ambient_coef'         : 0.4                  ,
        #  'light_shininess'            : 5.5                  ,
        #  'light_intensity'            : [0.6,0.6,0.6]        ,
        #  'light_specular_color'       : [1.0,1.0,1.0]        ,


        #light
        self.light_position       = glCore.light_position      
        self.light_color          = glCore.light_color         
        self.light_ambient_coef   = glCore.light_ambient_coef  
        self.light_shininess      = glCore.light_shininess     
        self.light_intensity      = glCore.light_intensity     
        self.light_specular_color = glCore.light_specular_color
        
        



        v_instances = """
        #version 330

        uniform mat4 model_mat;
        uniform mat4 view_mat;
        uniform mat4 proj_mat;

        in vec3 vert_coord;
        in vec3 vert_color;
        in vec3 vert_instance;
        in float vert_radius;

        vec3 vert_norm;

        out vec3 frag_coord;
        out vec3 frag_color;
        out vec3 frag_norm;

        void main(){
            mat4 modelview = view_mat * model_mat;
            vec3 offset_coord = vert_coord * vert_radius + vert_instance;
            gl_Position = proj_mat * modelview * vec4(offset_coord, 1.0);
            
            vert_norm = normalize(offset_coord - vert_instance);
            frag_coord = vec3(modelview * vec4(offset_coord, 1.0));
            frag_norm = mat3(transpose(inverse(model_mat))) * vert_norm;
            frag_color = vert_color;
        }
        """
        f_instances = """
        #version 330

        struct Light {
           vec3 position;
           //vec3 color;
           vec3 intensity;
           //vec3 specular_color;
           float ambient_coef;
           float shininess;
        };

        uniform Light my_light;

        uniform vec4 fog_color;
        uniform float fog_start;
        uniform float fog_end;



        in vec3 frag_coord;
        in vec3 frag_color;
        in vec3 frag_norm;

        out vec4 final_color;

        vec4 calculate_color(vec3 fnrm, vec3 fcrd, vec3 fcol){
            vec3 normal = normalize(fnrm);
            vec3 vert_to_light = normalize(my_light.position);
            vec3 vert_to_cam = normalize(fcrd);
            // Ambient Component
            vec3 ambient = my_light.ambient_coef * fcol * my_light.intensity;
            // Diffuse component
            float diffuse_coef = max(0.0, dot(normal, vert_to_light));
            vec3 diffuse = diffuse_coef * fcol * my_light.intensity;
            // Specular component
            float specular_coef = 0.0;
            if (diffuse_coef > 0.0)
                specular_coef = pow(max(0.0, dot(vert_to_cam, reflect(vert_to_light, normal))), my_light.shininess);
            vec3 specular = specular_coef * my_light.intensity;
            specular = specular * (vec3(1) - diffuse);
            vec4 out_color = vec4(ambient + diffuse + specular, 1.0);
            return out_color;
        }

        void main(){
            final_color = calculate_color(frag_norm, frag_coord, frag_color);
        
            float dist = abs(frag_coord.z);
            if(dist>=fog_start){
                float fog_factor = (fog_end-dist)/(fog_end-fog_start);
                final_color = mix(fog_color, final_color, fog_factor);
            }
            else{
                final_color = final_color;
                }
        }
        """


        sel_v_instances = """
        #version 330

        uniform mat4 model_mat;
        uniform mat4 view_mat;
        uniform mat4 proj_mat;

        in vec3 vert_coord;
        in vec3 vert_color;
        in vec3 vert_instance;
        in float vert_radius;

        vec3 vert_norm;

        out vec3 frag_coord;
        out vec3 frag_color;
        out vec3 frag_norm;

        void main(){
            mat4 modelview = view_mat * model_mat;
            vec3 offset_coord = vert_coord * vert_radius + vert_instance;
            gl_Position = proj_mat * modelview * vec4(offset_coord, 1.0);
            
            vert_norm = normalize(offset_coord - vert_instance);
            frag_coord = vec3(modelview * vec4(offset_coord, 1.0));
            frag_norm = mat3(transpose(inverse(model_mat))) * vert_norm;
            frag_color = vert_color;
        }
        """
        
        
        sel_f_instances = """
        #version 330

        struct Light {
           vec3 position;
           //vec3 color;
           vec3 intensity;
           //vec3 specular_color;
           float ambient_coef;
           float shininess;
        };

        uniform Light my_light;

        uniform vec4 fog_color;
        uniform float fog_start;
        uniform float fog_end;



        in vec3 frag_coord;
        in vec3 frag_color;
        in vec3 frag_norm;
        
        out vec4 final_color;

        vec4 calculate_color(vec3 fnrm, vec3 fcrd, vec3 fcol){
            vec3 normal = normalize(fnrm);
            vec3 vert_to_light = normalize(my_light.position);
            vec3 vert_to_cam = normalize(fcrd);
            // Ambient Component
            vec3 ambient = my_light.ambient_coef * fcol * my_light.intensity;
            // Diffuse component
            float diffuse_coef = max(0.0, dot(normal, vert_to_light));
            vec3 diffuse = diffuse_coef * fcol * my_light.intensity;
            // Specular component
            float specular_coef = 0.0;
            if (diffuse_coef > 0.0)
                specular_coef = pow(max(0.0, dot(vert_to_cam, reflect(vert_to_light, normal))), my_light.shininess);
            vec3 specular = specular_coef * my_light.intensity;
            specular = specular * (vec3(1) - diffuse);
            vec4 out_color = vec4(ambient + diffuse + specular, 1.0);
            return out_color;
        }

        void main(){
            //final_color = calculate_color(frag_norm, frag_coord, frag_color);
            final_color = vec4(frag_color, 1.0);
        
            //float dist = abs(frag_coord.z);
            //if(dist>=fog_start){
            //    float fog_factor = (fog_end-dist)/(fog_end-fog_start);
            //    final_color = mix(fog_color, final_color, fog_factor);
            //}
            //else{
            //    final_color = final_color;
            //    }
        }
        """




        self.gl_program_instances = self.load_shaders(v_instances, f_instances)
        self.gl_program_sel_instances = self.load_shaders(sel_v_instances, sel_f_instances)
        self.instances_vao = None
        self.insta_flag_test = None


    def load_lights(self, program):
        """ Function doc
        """
        light_pos = GL.glGetUniformLocation(program, "my_light.position")
        GL.glUniform3fv(light_pos, 1, self.light_position)
        #light_col = GL.glGetUniformLocation(program, "my_light.color")
        #GL.glUniform3fv(light_col, 1, self.light_color)
        amb_coef = GL.glGetUniformLocation(program, "my_light.ambient_coef")
        GL.glUniform1fv(amb_coef, 1, self.light_ambient_coef)
        shiny = GL.glGetUniformLocation(program, "my_light.shininess")
        GL.glUniform1fv(shiny, 1, self.light_shininess)
        intensity = GL.glGetUniformLocation(program, "my_light.intensity")
        GL.glUniform3fv(intensity, 1, self.light_intensity)
        #spec_col = GL.glGetUniformLocation(program, "my_light.specular_color")
        #GL.glUniform3fv(spec_col, 1, self.light_specular_color)
        return True

    def load_shaders(self, vertex, fragment, geometry=None):
        """ Here the shaders are loaded and compiled to an OpenGL program. By default
            the constructor shaders will be used, if you want to change the shaders
            use this function. The flag is used to create only one OpenGL program.
            
            Keyword arguments:
            vertex -- The vertex shader to be used
            fragment -- The fragment shader to be used
        """
        my_vertex_shader = self.create_shader(vertex, GL.GL_VERTEX_SHADER)
        my_fragment_shader = self.create_shader(fragment, GL.GL_FRAGMENT_SHADER)
        if geometry is not None:
            my_geometry_shader = self.create_shader(geometry, GL.GL_GEOMETRY_SHADER)
        program = GL.glCreateProgram()
        GL.glAttachShader(program, my_vertex_shader)
        GL.glAttachShader(program, my_fragment_shader)
        if geometry is not None:
            GL.glAttachShader(program, my_geometry_shader)
        GL.glLinkProgram(program)
        return program
    
    def create_shader(self, shader_prog, shader_type):
        """ Creates, links to a source, compiles and returns a shader.
            
            Keyword arguments:
            shader -- The shader text to use
            shader_type -- The OpenGL enum type of shader, it can be:
                           GL.GL_VERTEX_SHADER, GL.GL_GEOMETRY_SHADER or GL.GL_FRAGMENT_SHADER
            
            Returns:
            A shader object identifier or pops out an error
        """
        shader = GL.glCreateShader(shader_type)
        GL.glShaderSource(shader, shader_prog)
        GL.glCompileShader(shader)
        if GL.glGetShaderiv(shader, GL.GL_COMPILE_STATUS) != GL.GL_TRUE:
           #print("Error compiling the shader: ", shader_type)
            raise RuntimeError(GL.glGetShaderInfoLog(shader))
        return shader

        
    def _make_gl_vao_and_vbos (self, program, sel_program):
        ''' '''
        coords, indexes, colors = sphd.get_sphere([1,1,1], 1.0, [0, 1, 0], level="level_2")
        radii = np.ones(1, dtype=np.float32)
        instances = np.zeros(3,dtype=np.float32)
        #coords = coords.reshape(42,3)
        ##print(coords.shape)

        vao = GL.glGenVertexArrays(1)
        GL.glBindVertexArray(vao)

        ind_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indexes.nbytes, indexes, GL.GL_DYNAMIC_DRAW)

        coord_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.nbytes, coords, GL.GL_STATIC_DRAW)
        gl_coord = GL.glGetAttribLocation(program, "vert_coord")
        GL.glEnableVertexAttribArray(gl_coord)
        GL.glVertexAttribPointer(gl_coord, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))

        col_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_STATIC_DRAW)
        gl_colors = GL.glGetAttribLocation(program, "vert_color")
        GL.glEnableVertexAttribArray(gl_colors)
        GL.glVertexAttribPointer(gl_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))
        GL.glVertexAttribDivisor(gl_colors, 1)
        
        #sel_col_vbo = GL.glGenBuffers(1)
        #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, sel_col_vbo)
        #GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_STATIC_DRAW)
        #gl_colors = GL.glGetAttribLocation(sel_program, "vert_color")
        #GL.glEnableVertexAttribArray(gl_colors)
        #GL.glVertexAttribPointer(gl_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))
        #GL.glVertexAttribDivisor(gl_colors, 1)

        rad_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, rad_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, radii.nbytes, radii, GL.GL_STATIC_DRAW)
        gl_rads = GL.glGetAttribLocation(program, "vert_radius")
        GL.glEnableVertexAttribArray(gl_rads)
        GL.glVertexAttribPointer(gl_rads, 1, GL.GL_FLOAT, GL.GL_FALSE, radii.itemsize, ctypes.c_void_p(0))
        GL.glVertexAttribDivisor(gl_rads, 1)

        insta_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, insta_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, instances.nbytes, instances, GL.GL_STATIC_DRAW)
        gl_insta = GL.glGetAttribLocation(program, "vert_instance")
        GL.glEnableVertexAttribArray(gl_insta)
        GL.glVertexAttribPointer(gl_insta, 3, GL.GL_FLOAT, GL.GL_FALSE, 0, ctypes.c_void_p(0))
        GL.glVertexAttribDivisor(gl_insta, 1)

        GL.glBindVertexArray(0)
        GL.glDisableVertexAttribArray(gl_coord)
        GL.glDisableVertexAttribArray(gl_colors)
        GL.glDisableVertexAttribArray(gl_rads)
        GL.glDisableVertexAttribArray(gl_insta)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)

        return vao, (coord_vbo, col_vbo, rad_vbo, insta_vbo), int(len(indexes))
 
    
    def update_atomic_indexes (self, indexes = [] ):
        """ Function doc """
        self.indexes = indexes
        
        self.crd = []
        col     = []
        rads    = []
        sel_col = []
        
        
        for index in indexes:
            
            rads.append(self.visObj.vdw_dot_sizes[index]  )

            col.append(self.visObj.colors[index*3]  )
            col.append(self.visObj.colors[index*3+1])
            col.append(self.visObj.colors[index*3+2])
            
            sel_col.append(self.visObj.color_indexes[index*3]  )
            sel_col.append(self.visObj.color_indexes[index*3+1])
            sel_col.append(self.visObj.color_indexes[index*3+2])
        
        
        col       = np.array(col , dtype=np.float32)
        sel_col   = np.array(sel_col , dtype=np.float32)
        rads      = np.array(rads, dtype=np.float32)*self.scale
        
        for frame in self.visObj.frames:
            
            new_frame = []
            
            for index in indexes:
                
                new_frame.append(frame[index*3]  )
                new_frame.append(frame[index*3+1])
                new_frame.append(frame[index*3+2])
        
            new_frame = np.array(new_frame, dtype=np.float32)
            self.crd.append(new_frame)
    
        self.insta_col     = col
        self.insta_sel_col = sel_col
        self.insta_rads    = rads
        
        if self.instances_vao is None:
            self.instances_vao, self.instances_vbos, self.instances_elemns = self._make_gl_vao_and_vbos(self.gl_program_instances, self.gl_program_sel_instances )

        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.instances_vbos[1])
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.insta_col.nbytes, self.insta_col, GL.GL_STATIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.instances_vbos[2])
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.insta_rads.nbytes, self.insta_rads, GL.GL_STATIC_DRAW)
        
        
        
    def draw_representation (self):
        """ Function doc """
        
        if self.instances_vao is None:
            self.instances_vao, self.instances_vbos, self.instances_elemns = self._make_gl_vao_and_vbos(self.gl_program_instances, self.gl_program_sel_instances)
          
            self.update_atomic_indexes (indexes = self.indexes )
            
            '''
            self.insta_rads = self.visObj.vdw_dot_sizes*0.07
            self.insta_col  = np.array(self.visObj.colors, dtype=np.float32)
            self.insta_crd  = self.visObj.frames[0]
            '''
            self.insta_flag_test = True
            self.glCore.queue_draw()
        
        else:
            
            GL.glEnable(GL.GL_DEPTH_TEST)
            GL.glUseProgram(self.gl_program_instances)
            self.glCore.load_matrices(self.gl_program_instances, self.visObj.model_mat)
           
            self.load_lights(self.gl_program_instances)
            self.glCore.load_fog(self.gl_program_instances)
            
            GL.glBindVertexArray(self.instances_vao)
            if self.insta_flag_test:
                #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.instances_vbos[1])
                #GL.glBufferData(GL.GL_ARRAY_BUFFER, self.insta_col.nbytes, self.insta_col, GL.GL_STATIC_DRAW)
                #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.instances_vbos[2])
                #GL.glBufferData(GL.GL_ARRAY_BUFFER, self.insta_rads.nbytes, self.insta_rads, GL.GL_STATIC_DRAW)
                #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.instances_vbos[3])
                #GL.glBufferData(GL.GL_ARRAY_BUFFER, self.insta_crd.nbytes, self.insta_crd, GL.GL_STATIC_DRAW)
                self.insta_flag_test = False
            
            
            self.insta_crd = self.crd[self.glCore._safe_frame_exchange(visObj = self.visObj, return_frame = False)]#self.glCore._safe_frame_exchange(self.visObj)


            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.instances_vbos[1])
            GL.glBufferData(GL.GL_ARRAY_BUFFER, self.insta_col.nbytes, self.insta_col, GL.GL_STATIC_DRAW)

            
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.instances_vbos[3])
            GL.glBufferData(GL.GL_ARRAY_BUFFER, self.insta_crd.nbytes, self.insta_crd, GL.GL_STATIC_DRAW)
            
            #print (frame)
            GL.glDrawElementsInstanced(GL.GL_TRIANGLES, self.instances_elemns, GL.GL_UNSIGNED_INT, None, self.insta_crd.shape[0])
            GL.glBindVertexArray(0)
            GL.glUseProgram(0)        
            GL.glDisable(GL.GL_DEPTH_TEST)

    def draw_background_sel_representation  (self, line_width_factor = 5):
        """ Function doc """
        
        if self.instances_vao is None:
            self.instances_vao, self.instances_vbos, self.instances_elemns = self._make_gl_vao_and_vbos(self.gl_program_instances)
          
            self.update_atomic_indexes (indexes = self.indexes )
            
            '''
            self.insta_rads = self.visObj.vdw_dot_sizes*0.07
            self.insta_col  = np.array(self.visObj.colors, dtype=np.float32)
            self.insta_crd  = self.visObj.frames[0]
            '''
            self.insta_flag_test = True
            self.glCore.queue_draw()
        
        else:
            
            GL.glEnable(GL.GL_DEPTH_TEST)
            GL.glUseProgram(self.gl_program_sel_instances)
            self.glCore.load_matrices(self.gl_program_sel_instances, self.visObj.model_mat)
           
            self.load_lights(self.gl_program_sel_instances)
            self.glCore.load_fog(self.gl_program_sel_instances)
            
            GL.glBindVertexArray(self.instances_vao)
            if self.insta_flag_test:
                #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.instances_vbos[1])
                #GL.glBufferData(GL.GL_ARRAY_BUFFER, self.insta_col.nbytes, self.insta_col, GL.GL_STATIC_DRAW)
                #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.instances_vbos[2])
                #GL.glBufferData(GL.GL_ARRAY_BUFFER, self.insta_rads.nbytes, self.insta_rads, GL.GL_STATIC_DRAW)
                #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.instances_vbos[3])
                #GL.glBufferData(GL.GL_ARRAY_BUFFER, self.insta_crd.nbytes, self.insta_crd, GL.GL_STATIC_DRAW)
                self.insta_flag_test = False
            
            
            self.insta_crd = self.crd[self.glCore._safe_frame_exchange(visObj = self.visObj, return_frame = False)]#self.glCore._safe_frame_exchange(self.visObj)
           
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.instances_vbos[1])
            GL.glBufferData(GL.GL_ARRAY_BUFFER, self.insta_sel_col.nbytes, self.insta_sel_col, GL.GL_STATIC_DRAW)
            
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.instances_vbos[3])
            GL.glBufferData(GL.GL_ARRAY_BUFFER, self.insta_crd.nbytes, self.insta_crd, GL.GL_STATIC_DRAW)
            
            #print (frame)
            GL.glDrawElementsInstanced(GL.GL_TRIANGLES, self.instances_elemns, GL.GL_UNSIGNED_INT, None, self.insta_crd.shape[0])
            GL.glBindVertexArray(0)
            GL.glUseProgram(0)        
            GL.glDisable(GL.GL_DEPTH_TEST)













class GlumpyRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'glumpy', active = True, _type = 'mol', visObj = None, glCore = None, scale=1.0):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type
        self.visObj             = visObj
        self.glCore             = glCore
        self.scale              = scale
        self.vbos_list =[]

        # representation    
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None

        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None


    def _make_gl_vao_and_vbos (self, indexes = None):
        """ Function doc """
        
        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        
        coords  = self.visObj.frames[0]
        colors  = self.visObj.colors
        radii   = [self.scale] * len(self.visObj.frames[0])
        radii   = np.array(radii)
        ##print ('radii', radii)

        dot_qtty  = int(len(coords)/3)
        indexes = np.arange(dot_qtty, dtype=np.uint32)

        self._make_gl_representation_vao_and_vbos (indexes    = indexes,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None  ,
                                                   )
        colors_idx = self.visObj.color_indexes
        self._make_gl_sel_representation_vao_and_vbos (indexes    = indexes    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None      ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs()
        self._enable_anti_alis_to_lines()
        GL.glUseProgram(self.shader_program)
        #GL.glEnable(GL.GL_VERTEX_PROGRAM_POINT_SIZE)
        height = self.visObj.vm_session.glwidget.vm_widget.height
        dist_cam_zrp = self.visObj.vm_session.glwidget.vm_widget.dist_cam_zrp
        
        #GL.glPointSize((50*height/(abs(dist_cam_zrp)))**0.5)
        #GL.glPointSize(55)
        #print('passei aqui')
        
        xyz_coords = self.glCore.glcamera.get_modelview_position(self.visObj.model_mat)
        u_campos = GL.glGetUniformLocation(self.shader_program, 'u_campos')
        GL.glUniform3fv(u_campos, 1, xyz_coords)
        
        #u_depth = GL.glGetUniformLocation(self.shader_program, 'u_depth')
        #GL.glUniform1fv(u_depth, 1, (self.glCore.glcamera.z_near - self.glCore.glcamera.z_far))
        
        
        self.glCore.load_lights  (self.shader_program )

        
        
        
        
        self.glCore.load_matrices(self.shader_program, self.visObj.model_mat)
        self.glCore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        if self.glCore.modified_view:
            pass
        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
            GL.glDrawElements(GL.GL_POINTS, int(len(self.visObj.atoms)), GL.GL_UNSIGNED_INT, None)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glUseProgram(self.sel_shader_program)
        self.glCore.load_matrices(self.sel_shader_program, self.visObj.model_mat)
        GL.glBindVertexArray(self.sel_vao)
        if self.glCore.modified_view:
            pass
        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea
            '''
            #print(self.name,'draw_background_sel_representation')
            self._set_coordinates_to_buffer (coord_vbo = False, sel_coord_vbo = True)
            GL.glDrawElements(GL.GL_POINTS, int(len(self.visObj.atoms)), GL.GL_UNSIGNED_INT, None)









class CartoonRepresentation (Representation):
    def __init__ (self, name = 'cartoon', active = True, _type = 'mol', visObj = None, glCore = None, indexes = []):
        self.name               = name
        self.active             = active
        self.type               = _type

        self.visObj             = visObj
        self.glCore             = glCore
        
        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.norm_vbo       = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        
        
        coords, normals, indexes, colors = cartoon.cartoon(visObj, spline_detail=5)
        
        coords = coords.flatten()
        normals = normals.flatten()
        colors = colors.flatten()
        
        
        self.coords2 = coords
        self.colors2 = colors
        self.normals2 = normals
        self.indexes2 = indexes
        self.vbos_list =[]

    def _make_gl_vao_and_vbos (self, indexes = None):
        """ Function doc """
        #if indexes is not None:
        #    pass
        #else:
        
        #dot_qtty  = int(len(self.visObj.frames[0])/3)
        #indexes = []
        #for i in range(dot_qtty):
        #    indexes.append(i)
        

        self.shader_program     = self.glCore.shader_programs[self.name]
        #self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        

        '''
        coords  = np.array(self.coords2, dtype=np.float32)
        colors  = np.array(self.colors2, dtype=np.float32)
        normals = np.array(self.normals2, dtype=np.float32)
        indexes = np.array(self.indexes2, dtype=np.uint32)
        '''
        
        
        coords  = self.coords2 
        colors  = self.colors2 
        normals = self.normals2
        indexes = self.indexes2
        
        #print ('len(coords),len(colors), len(normals),len(indexes)', len(coords),len(colors), len(normals),len(indexes)  )

        self._make_gl_representation_vao_and_vbos (indexes    = indexes,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   normals    = normals
                                                   )
        
        
        
        self.ind_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indexes.itemsize*len(indexes), indexes, GL.GL_DYNAMIC_DRAW)
        
        #self.coord_vbo = GL.glGenBuffers(1)
        #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.coord_vbo)
        ##GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.itemsize*len(coords), coords, GL.GL_STATIC_DRAW)
        #GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.nbytes, coords, GL.GL_STATIC_DRAW)
        #gl_coord = GL.glGetAttribLocation(self.shader_program, 'vert_coord')
        #GL.glEnableVertexAttribArray(gl_coord)
        #GL.glVertexAttribPointer(gl_coord, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
        
        
        self.col_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.itemsize*len(colors), colors, GL.GL_STATIC_DRAW)
        gl_color = GL.glGetAttribLocation(self.shader_program, 'vert_color')
        GL.glEnableVertexAttribArray(gl_color)
        GL.glVertexAttribPointer(gl_color, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))

        self.norm_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.norm_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, normals.itemsize*len(normals), normals, GL.GL_STATIC_DRAW)
        gl_norm = GL.glGetAttribLocation(self.shader_program, 'vert_norm')
        GL.glEnableVertexAttribArray(gl_norm)
        GL.glVertexAttribPointer(gl_norm, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*normals.itemsize, ctypes.c_void_p(0))
        
        
        
        
        
        #self.centr_vbo = GL.glGenBuffers(1)
        #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coords)
        #GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.itemsize*len(coords), coords, GL.GL_STATIC_DRAW)
        #gl_center = GL.glGetAttribLocation(self.shader_program , 'vert_centr')
        #GL.glEnableVertexAttribArray(gl_center)
        #GL.glVertexAttribPointer(gl_center, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
        
        
        
        colors_idx = self.visObj.color_indexes
        self.sel_vao = True
        '''
        self._make_gl_sel_representation_vao_and_vbos (indexes    = indexes    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )
        '''
    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        #self._enable_anti_alis_to_lines()
        
        
        
        
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glDisable(GL.GL_CULL_FACE)
        #GL.glCullFace(GL.GL_BACK)
        view = self.glCore.glcamera.view_matrix
        
        GL.glUseProgram(self.shader_program )
        
        #print (self.visObj.model_mat,view)
        
        m_normal = np.array(np.matrix(np.dot(view, self.visObj.model_mat)).I.T)
        
        self.glCore.load_matrices(self.shader_program , self.visObj.model_mat)
        self.glCore.load_lights  (self.shader_program )
        self.glCore.load_fog     (self.shader_program )
        GL.glBindVertexArray(self.vao)
        
        
        
        
        
        
        
        
        
        '''
        #print ('DotsRepresentation')
        height = self.visObj.vm_session.glwidget.vm_widget.height
        
        GL.glUseProgram(self.shader_program)
        #1*self.height dot_size
        #GL.glLineWidth(40/abs(self.glCore.dist_cam_zrp))
        GL.glPointSize(0.1*height/abs(self.glCore.dist_cam_zrp)) # dot size not included yet
        self.glCore.load_matrices(self.shader_program, self.visObj.model_mat)
        self.glCore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        '''
        if self.glCore.modified_view:
            pass
        
        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            # self._set_coordinates_to_buffer(coord_vbo = True, sel_coord_vbo = False)
            #GL.glDrawElements(GL.GL_POINTS, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
            #GL.glDrawElements(GL.GL_LINE_LOOP, int(len(self.coords2)), GL.GL_UNSIGNED_INT, None)
            #GL.glDrawElements(GL.GL_LINE_STRIP, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
            
            #print('int(len(self.indexes2))', int(len(self.indexes2)))
            GL.glDrawElements(GL.GL_TRIANGLES, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
            #GL.glDrawElements(GL.GL_TRIANGLES, 54060, GL.GL_UNSIGNED_INT, None)
        
        #GL.glBindVertexArray(0)
        #GL.glLineWidth(1)
        #GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        pass






        
class SurfaceRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'surface', active = True, _type = 'mol', visObj = None, glCore = None, indexes = []):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type

        self.visObj             = visObj
        self.glCore             = glCore
        
        self.vbos_list =[]


        
        
        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.norm_vbo       = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        self.read_surface_data()
    
    
    ##### sub 2 vev3 vectors
    def sub_vec3(self, a, b):
        c = [ a[0] - b[0],
              a[1] - b[1],
              a[2] - b[2] ]

        return c

    ## add 2 vectors and take the avg
    ## if a vector is still 0 we just take b
    def avg_add_vec3(self, a, b):
        if a[0] == 0.0 and a[1] == 0.0 and a[2] == 0.0 :
            return b

        c = [ (a[0] + b[0]) * 0.5 ,
              (a[1] + b[1]) * 0.5 ,
              (a[2] + b[2]) * 0.5 ]

        return c    

    ## make the cross product of 2 vectors
    def cross_vec3(self, a, b):
        c = [a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0]]

        return c
    #############################################
        
    
    
    def read_surface_data(self):
        """ Function doc """
        #from random import random 
        #
        #[verts, tris, verts_gpu, tris_gpu] = edtsurf.calc_surface('/home/fernando/programs/EasyHybrid3/Coords/pdbs/1bx4_H.pdb')
        #self.coords2  = verts_gpu
        #self.indexes2 = tris_gpu
        #self.colors2  = []
        #
        #
        #size = len( self.coords2 )
        #for i in range(size):
        #    self.colors2.append(float(i/size) + random())
        
        rawdata = open('../EasyHybrid3/Coords/pdbs/1bx4.ply', 'r')
        lines  = rawdata.readlines()
        
        self.coords2 = []
        self.colors2 = []
        self.normals2 = []
        self.indexes2 = []
        avg_normals_indexes = []
        
        
        for line in lines:
            line2 = line.split()
            
            if len(line2) == 6:
                #print (line2)
                self.coords2.append(float(line2[0]))
                self.coords2.append(float(line2[1]))
                self.coords2.append(float(line2[2]))
                                                  
                self.colors2.append(float(line2[3])/255)
                self.colors2.append(float(line2[4])/255)
                self.colors2.append(float(line2[5])/255)
                
                self.normals2.append(float(line2[0]))
                self.normals2.append(float(line2[1]))
                self.normals2.append(float(line2[2]))                
                avg_normals_indexes.append( ( 0.0 , 0.0 , 0.0 ) )  ### NEW !!! 

            if len(line2) == 7:
                
                self.indexes2.append(int(line2[1]))
                self.indexes2.append(int(line2[2]))
                self.indexes2.append(int(line2[3]))
                
        
        ## calculate normals and interpolate them (thanks a lot Kai)
        for i in range( 0 , len(self.indexes2) , 3 ):

            index_1 = self.indexes2[i] * 3;
            index_2 = self.indexes2[i+1] * 3;
            index_3 = self.indexes2[i+2] * 3;
            vertex_1 = ( self.coords2[index_1] , self.coords2[index_1+1] , self.coords2[index_1+2] )
            vertex_2 = ( self.coords2[index_2] , self.coords2[index_2+1] , self.coords2[index_2+2] )
            vertex_3 = ( self.coords2[index_3] , self.coords2[index_3+1] , self.coords2[index_3+2] )

            vec_p0_p1 = self.sub_vec3( vertex_2 , vertex_1 )
            vec_p0_p2 = self.sub_vec3( vertex_3 , vertex_1 )
            norm_vec  = self.cross_vec3( vec_p0_p1, vec_p0_p2 )

            vert_index_1 = self.indexes2[i] ;
            vert_index_2 = self.indexes2[i+1] ;
            vert_index_3 = self.indexes2[i+2] ;
            
            avg_normals_indexes[vert_index_1] = self.avg_add_vec3( avg_normals_indexes[vert_index_1] , norm_vec )
            avg_normals_indexes[vert_index_2] = self.avg_add_vec3( avg_normals_indexes[vert_index_2] , norm_vec )
            avg_normals_indexes[vert_index_3] = self.avg_add_vec3( avg_normals_indexes[vert_index_3] , norm_vec )


        ## set all new interpolated normals   
        for i in range( 0 , len(self.indexes2) , 1 ):
            index_1 = self.indexes2[i] * 3;

            self.normals2[index_1]   = avg_normals_indexes[self.indexes2[i]][0]
            self.normals2[index_1+1] = avg_normals_indexes[self.indexes2[i]][1]
            self.normals2[index_1+2] = avg_normals_indexes[self.indexes2[i]][2]





               
                

    def _make_gl_vao_and_vbos (self, indexes = None):
        """ Function doc """
        #if indexes is not None:
        #    pass
        #else:
        
        #dot_qtty  = int(len(self.visObj.frames[0])/3)
        #indexes = []
        #for i in range(dot_qtty):
        #    indexes.append(i)
        

        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        
        #indexes = np.array(self.visObj.index_bonds, dtype=np.uint32)
        #indexes = np.array(self.visObj.idex, dtype=np.uint32)

        coords  = np.array(self.coords2, dtype=np.float32)
        colors  = np.array(self.colors2, dtype=np.float32)
        normals = np.array(self.normals2, dtype=np.float32)
        #indexes = range(0, len(self.coords2))     
        #indexes = np.array(indexes, dtype=np.uint32)
        indexes = np.array(self.indexes2, dtype=np.uint32)
        #print (indexes)


        self._make_gl_representation_vao_and_vbos (indexes    = indexes,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   normals    = normals
                                                   )
        
        #self.centr_vbo = GL.glGenBuffers(1)
        #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coords)
        #GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.itemsize*len(coords), coords, GL.GL_STATIC_DRAW)
        #gl_center = GL.glGetAttribLocation(self.shader_program , 'vert_centr')
        #GL.glEnableVertexAttribArray(gl_center)
        #GL.glVertexAttribPointer(gl_center, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
        
        
        
        colors_idx = self.visObj.color_indexes
        self._make_gl_sel_representation_vao_and_vbos (indexes    = indexes    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        #self._enable_anti_alis_to_lines()
        
        
        
        
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_CULL_FACE)
        GL.glCullFace(GL.GL_BACK)
        view = self.glCore.glcamera.view_matrix
        
        GL.glUseProgram(self.shader_program )
        
        #print (self.visObj.model_mat,view)
        
        m_normal = np.array(np.matrix(np.dot(view, self.visObj.model_mat)).I.T)
        
        self.glCore.load_matrices(self.shader_program , self.visObj.model_mat)
        self.glCore.load_lights  (self.shader_program )
        self.glCore.load_fog     (self.shader_program )
        GL.glBindVertexArray(self.vao)
        
        
        
        
        
        
        
        
        
        '''
        #print ('DotsRepresentation')
        height = self.visObj.vm_session.glwidget.vm_widget.height
        
        GL.glUseProgram(self.shader_program)
        #1*self.height dot_size
        #GL.glLineWidth(40/abs(self.glCore.dist_cam_zrp))
        GL.glPointSize(0.1*height/abs(self.glCore.dist_cam_zrp)) # dot size not included yet
        self.glCore.load_matrices(self.shader_program, self.visObj.model_mat)
        self.glCore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        '''
        if self.glCore.modified_view:
            pass
        
        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            # self._set_coordinates_to_buffer(coord_vbo = True, sel_coord_vbo = False)
            #GL.glDrawElements(GL.GL_POINTS, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
            #GL.glDrawElements(GL.GL_LINE_LOOP, int(len(self.coords2)), GL.GL_UNSIGNED_INT, None)
            #GL.glDrawElements(GL.GL_LINE_STRIP, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
            GL.glDrawElements(GL.GL_TRIANGLES, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
        
        #GL.glBindVertexArray(0)
        #GL.glLineWidth(1)
        #GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        pass





class WiresRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'wires', active = True, _type = 'mol', visObj = None, glCore = None, indexes = []):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type
        self.visObj             = visObj
        self.glCore             = glCore
        self.vbos_list =[]

        # representation    
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
        
        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None

        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        self.read_surface_data()
    
    def read_surface_data(self):
        """ Function doc """
        rawdata = open('../EasyHybrid3/Coords/pdbs/1bx4.ply', 'r')
        lines  = rawdata.readlines()
        
        self.coords2 = []
        self.colors2 = []
        self.indexes2 = []
        
        for line in lines:
            line2 = line.split()
            if len(line2) == 6:
                self.coords2.append(float(line2[0]))
                self.coords2.append(float(line2[1]))
                self.coords2.append(float(line2[2]))
                self.colors2.append(float(line2[3])/255)
                self.colors2.append(float(line2[4])/255)
                self.colors2.append(float(line2[5])/255)
            if len(line2) == 7:
                self.indexes2.append(int(line2[1]))
                self.indexes2.append(int(line2[2]))
                self.indexes2.append(int(line2[3]))

    def _make_gl_vao_and_vbos (self, indexes = None):
        """ Function doc """
        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        coords  = np.array(self.coords2, dtype=np.float32)
        colors  = np.zeros(len(self.colors2))
        indexes = np.array(self.indexes2, dtype=np.uint32)
        self._make_gl_representation_vao_and_vbos (indexes    = indexes,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   )
        colors_idx = self.visObj.color_indexes
        self._make_gl_sel_representation_vao_and_vbos (indexes    = indexes    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        pass
        #GL.glEnable(GL.GL_DEPTH_TEST)
        #GL.glEnable(GL.GL_CULL_FACE)
        #GL.glCullFace(GL.GL_BACK)
        #
        ##LineWidth = (80/abs(self.glCore.dist_cam_zrp)/2)**0.5  #40/abs(self.glCore.dist_cam_zrp)
        ##GL.glLineWidth(2)
        #
        #
        #view = self.glCore.glcamera.view_matrix
        #GL.glUseProgram(self.shader_program )
        #m_normal = np.array(np.matrix(np.dot(view, self.visObj.model_mat)).I.T)
        #self.glCore.load_matrices(self.shader_program , self.visObj.model_mat)
        ##self.glCore.load_lights  (self.shader_program )
        #self.glCore.load_fog     (self.shader_program )
        #GL.glBindVertexArray(self.vao)
        #if self.glCore.modified_view:
        #    pass
        #
        #else:
        #    '''
        #    This function checks if the number of the called frame will not exceed 
        #    the limit of frames that each object has. Allowing two objects with 
        #    different trajectory sizes to be manipulated at the same time within the 
        #    glArea'''
        #    # self._set_coordinates_to_buffer(coord_vbo = True, sel_coord_vbo = False)
        #    GL.glDrawElements(GL.GL_TRIANGLES, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
        #GL.glDisable(GL.GL_DEPTH_TEST)
        
    def draw_background_sel_representation  (self):
        """ Function doc """
        pass







class LabelRepresentation:
    """ Class doc """
    
    def __init__ (self, name = 'labels', active = True, _type = 'mol', visObj = None, glCore = None, indexes = []):
        """ Class initialiser """
        self.visObj = visObj
        self.name   = name
        self.active = True
        self.glCore = glCore
        self.vbos_list =[]

        self.chars     = 0 
        #self._check_VAO_and_VBOs()
        
    def _check_VAO_and_VBOs (self, indexes = None):
        """ Function doc """
        if self.visObj.vm_font.vao is None:
            self.visObj.vm_font.make_freetype_font()
            self.visObj.vm_font.make_freetype_texture(self.glCore.freetype_program)
        
        if self.chars == 0:
           #print('self._build_buffer()')
            self._build_buffer()
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.visObj.vm_font.vbos[0])
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.xyz_pos.itemsize*len(self.xyz_pos), self.xyz_pos, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.visObj.vm_font.vbos[1])
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.uv_coords.itemsize*len(self.uv_coords), self.uv_coords, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)


    def _build_buffer (self, indexes = None):
        self.chars     = 0
        self.xyz_pos   = []
        self.uv_coords = []
        for atom in self.visObj.atoms:
            
            texto = atom.name
            point = np.array(atom.coords (self.glCore.frame),np.float32)
            point = np.array((point[0],point[1],point[2],1),np.float32)
            point = np.dot(point, self.visObj.model_mat)

            GL.glBindTexture(GL.GL_TEXTURE_2D, self.visObj.vm_font.texture_id)
            for i,c in enumerate(texto):
                self.chars += 1
                c_id = ord(c)
                x = c_id%16
                y = c_id//16-2
                self.xyz_pos.append(point[0]+i*self.visObj.vm_font.char_width)
                self.xyz_pos.append(point[1])
                self.xyz_pos.append(point[2])

                self.uv_coords.append(x*self.visObj.vm_font.text_u)
                self.uv_coords.append(y*self.visObj.vm_font.text_v)
                self.uv_coords.append((x+1)*self.visObj.vm_font.text_u)
                self.uv_coords.append((y+1)*self.visObj.vm_font.text_v)
            #print(texto)
        #print('xyz_pos  ',len(self.xyz_pos))
        #print('uv_coords',len(self.uv_coords))
        #print('atoms    ',len(self.visObj.atoms))
        #print('chars    ',self.chars)
        
        self.xyz_pos   = np.array(self.xyz_pos  , np.float32)
        self.uv_coords = np.array(self.uv_coords, np.float32)
        


    
    
    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs()
        
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        GL.glUseProgram(self.glCore.freetype_program)
        
        self.visObj.vm_font.load_matrices(self.glCore.freetype_program, self.glCore.glcamera.view_matrix, self.glCore.glcamera.projection_matrix)
        self.visObj.vm_font.load_font_params(self.glCore.freetype_program)
        
        GL.glBindVertexArray(self.visObj.vm_font.vao)
        GL.glDrawArrays(GL.GL_POINTS, 0, self.chars)
        GL.glDisable(GL.GL_BLEND)
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)

        

    def draw_background_sel_representation  (self):
        """ Function doc """
        pass



