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
import gc
import os


gi.require_version('Gdk', '3.0')
from gi.repository import Gtk, Gdk
from matplotlib.backends.backend_gtk3agg import FigureCanvas  # or gtk3cairo.
from matplotlib.figure import Figure
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors





VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')


class PotentialEnergyAnalysisWindow():
    
    def OpenWindow (self, vobject = None):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/PES_analysis_window.glade'))
            self.builder.connect_signals(self)
            self.vobject = vobject
            self.window = self.builder.get_object('window')
            self.window.set_title('Analysis Window')
            self.window.set_keep_above(True)            
            self.window.set_default_size(700, 450)
            
            self.hbox           = self.builder.get_object('hbox_matplot_figures')
            self.scale_traj     = self.builder.get_object('scale_trajectory_from_PES')
            self.adjustment     = Gtk.Adjustment(value         = 0,
                                                 lower         = 0,
                                                 upper         = 100,
                                                 step_increment= 1,
                                                 page_increment= 1,
                                                 page_size     = 1)



            self.scale_traj.set_adjustment ( self.adjustment)
            self.scale_traj.set_digits(0)
            
            
            
            
            
            
            
            self.fig = Figure(figsize=(6, 4),constrained_layout=True)
            self.canvas = FigureCanvas(self.fig)  # a Gtk.DrawingArea
            self.canvas.mpl_connect('button_press_event', self.onpick2)

            self.hbox.pack_start(self.canvas, True, True, 0)
            self.ax = self.fig.add_subplot()
            
            
            
            
            sys_selected = 0
            #'''------------------------------------------------------------------------------------
            self.coordinates_combobox = self.builder.get_object('coordinates_combobox')
            renderer_text = Gtk.CellRendererText()
            self.coordinates_combobox.add_attribute(renderer_text, "text", 0)
            

            self.vobject_liststore = Gtk.ListStore(str, int)
            names = [ ]
            for key , system in self.main.pDynamo_session.systems.items():
                
                for vobject_id in  system['logfile_data'].keys():
                    _vobject = self.main.vm_session.vismol_objects_dic[vobject_id]
                    print(['_vobject:', _vobject.name, int(key)])
                    self.vobject_liststore.append([_vobject.name, int(key)])
                    

            self.coordinates_combobox.set_model(self.vobject_liststore)

            if sys_selected:
                self.coordinates_combobox.set_active(sys_selected)
            else:
                self.coordinates_combobox.set_active(0)
            #------------------------------------------------------------------------------------'''



            
            if self.vobject:
                data =self.main.pDynamo_session.systems[self.vobject.easyhybrid_system_id]['logfile_data'][self.vobject.index][0][1]
                #data = self.vobject.trajectory2D_data
            else:
                data = None#parse_2D_scan_logfile (logfile = '/home/fernando/programs/VisMol/examples/tobachega/ScanTraj.log')
            

            
            
            self.Z = data['Z']


            self.X = range(len(self.Z[0]))
            self.Y = range(len(self.Z))
            print(self.X)
            print(self.Y)
            print(self.Z)
            #self.pcm = self.ax.pcolormesh(self.X, self.Y, self.Z, cmap='jet', vmin=0, shading='gouraud')
            self.pcm = self.ax.pcolormesh(self.X, self.Y, self.Z, cmap='jet', vmin=0)#, shading='gouraud')
            am = self.ax.contour(self.X, self.Y, self.Z,lspacing = 10, colors='k')
            
            self.fig.colorbar(self.pcm,  ax=self.ax)#, extend='both')
            

            '''----------------------------------------------------------------------------------------------------'''
            
            fig = Figure(figsize=(5, 4), dpi=100)
            self.ax2 = fig.add_subplot( )
            self.ax3 = fig.add_subplot( )
            
            self.line2, = self.ax2.plot([], [], '-o')
            secax = self.ax2.secondary_xaxis('top', functions=(None))

            self.canvas2 = FigureCanvas(fig)  # a Gtk.DrawingArea
            self.hbox.pack_start(self.canvas2, True, True, 0)        
            
 
            self.window.show_all()
            self.Visible  = True
    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
    
    def onpick2(self, event):
        ''' Function doc '''
        self.ax.cla()
        #self.ax.pcolormesh(self.X, self.Y, self.Z, cmap='jet', vmin=0, shading='gouraud')
        #self.ax.pcolormesh(self.X, self.Y, self.Z, cmap='jet', vmin=0, shading='gouraud')
        self.pcm = self.ax.pcolormesh(self.X, self.Y, self.Z, cmap='jet', vmin=0)#, shading='gouraud')
        am = self.ax.contour(self.X, self.Y, self.Z,lspacing = 10, colors='k')
        if event.xdata is None or event.ydata is None:
            self.xdata  = []
            self.ydata  = []
            self.indexes = []
            self.zdata   = []
            self.ax2.cla()
            
            self.xy_traj = []
            self.scale_traj_new_definitions()

        else:
            
            print('you pressed', event.button)
            print('you pressed', event.xdata, event.ydata)
            x, y = int(event.xdata), int(event.ydata)
            print (x, y)
            print('you pressed', x, y, self.Z[y][x])
            print('you pressed', event)
            #self.xdata.append(event.xdata)
            #self.ydata.append(event.ydata)
            self.xdata.append(x+0.5)
            self.ydata.append(y+0.5)
            self.xy_traj.append([x,y])
            
            self.scale_traj_new_definitions()
            
            self.zdata.append(self.Z[y][x])
            #self.line, = self.ax.plot(event.xdata, event.ydata, '-o')  # plot the first row
        
        self.line, = self.ax.plot(self.xdata, self.ydata, '-ok')  # plot the first row
        self.line2, = self.ax2.plot(range(0, len(self.zdata)), self.zdata, '-ob')  # plot the first row
        
        self.canvas.draw()
        self.canvas2.draw()
        print()  

    def __init__(self, main = None ):
        """ Class initialiser """
        self.main     = main
        #self.p_session           = self.easyhybrid_main.pDynamo_session
        #self.vm_session          = main.vm_session
        self.Visible             =  False        
        
        self.vobject_liststore   = Gtk.ListStore(str, int)
        self.data_liststore      = Gtk.ListStore(str, int)
        #self.residue_liststore   = Gtk.ListStore(str, str, str)
        
        
        self.opt_methods = { 0 : 'ConjugatedGradient',
                             1 : 'SteepestDescent'   ,
                             2 : 'LFBGS'             ,
                             3 : 'QuasiNewton'       ,
                             4 : 'FIRE'              }
        
        self.xdata = []
        self.ydata = []
        self.zdata = []
        self.xy_traj = []

        self.vobject = None

    def scale_traj_new_definitions(self):
        #self.scale_traj
        self.scale_traj.set_range(0, len(self.xy_traj))
        #self.scale_traj.set_value(self.vm_session.get_frame())

    def on_scaler_frame_value_changed (self, hscale, text= None,  data=None):
        """ Function doc """
        value = self.scale_traj.get_value()
        pos   = self.scale_traj.get_value_pos ()
        
        self.ax3.cla()
        
        #print(self.xy_traj[int(value)])
        xy = self.xy_traj[int(value)]
        print(xy, self.zdata[int(value)])
        
        
        self.ax2.plot(range(0, len(self.zdata)), self.zdata, '-ob')
        self.ax3.plot( [int(value)], [self.zdata[int(value)]], '-or')
        
        if self.vobject:
            frame = self.vobject.trajectory2D_xy_indexes[(xy[0], xy[1])]
            self.main.vm_session.set_frame(int(frame)) 
        self.canvas2.draw()
        #self.scale_traj.set_value(int(value))
        
    def change_cb_coordType1 (self, combo_box):
        """ Function doc """
        
        _type = self.combobox_reaction_coord1.get_active()        
        
        if _type == 0:
            self.builder.get_object('label_atom3_coord1').hide()
            self.builder.get_object('entry_atom3_index_coord1').hide()
            self.builder.get_object('label_name3_coord1').hide()
            self.builder.get_object('entry_atom3_name_coord1').hide()
            
            self.builder.get_object('label_atom4_coord1').hide()
            self.builder.get_object('entry_atom4_index_coord1').hide()
            self.builder.get_object('label_name4_coord1').hide()
            self.builder.get_object('entry_atom4_name_coord1').hide()
            self.builder.get_object('mass_restraints1').set_sensitive(False)

        if _type == 1:
            self.builder.get_object('label_atom3_coord1').show()
            self.builder.get_object('entry_atom3_index_coord1').show()
            self.builder.get_object('label_name3_coord1').show()
            self.builder.get_object('entry_atom3_name_coord1').show()
            
            self.builder.get_object('label_atom4_coord1').hide()
            self.builder.get_object('entry_atom4_index_coord1').hide()
            self.builder.get_object('label_name4_coord1').hide()
            self.builder.get_object('entry_atom4_name_coord1').hide()
            self.builder.get_object('mass_restraints1').set_sensitive(True)

        if _type == 2:
            self.builder.get_object('label_atom3_coord1').show()
            self.builder.get_object('entry_atom3_index_coord1').show()
            self.builder.get_object('label_name3_coord1').show()
            self.builder.get_object('entry_atom3_name_coord1').show()
            
            self.builder.get_object('label_atom4_coord1').show()
            self.builder.get_object('entry_atom4_index_coord1').show()
            self.builder.get_object('label_name4_coord1').show()
            self.builder.get_object('entry_atom4_name_coord1').show()
            self.builder.get_object('mass_restraints1').set_sensitive(False)

        try:
            self.refresh_dmininum ( coord1 = True)
        except:
            print(texto_d1)
            print(texto_d2d1)
                    
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
    def run_scan(self,button):
        '''
        Get infotmation and run the simulation
        '''         
        parameters = {"simulation_type":"Relaxed_Surface_Scan",
                      "ndim":1                                ,
                      "ATOMS_RC1":None                        ,
                      "ATOMS_RC2":None                        ,
                      "nsteps_RC1":0                          ,
                      "nsteps_RC2":0                          ,
                      "force_constant_1":4000.0               ,
                      "force_constant_2":4000.0               ,
                      "maxIterations":1000                    ,
                      "dincre_RC1":0.1                        ,
                      "dincre_RC2":0.1                        ,
                      "dminimum_RC1":0.0                      ,
                      "dminimum_RC2":0.0                      ,
                      "sigma_pk1pk3_rc1":1.0                  ,
                      "sigma_pk3pk1_rc1":-1.0                 ,
                      "sigma_pk1pk3_rc2":1.0                  ,
                      "sigma_pk3pk1_rc2":-1.0                 ,
                      "rc_type_1":"Distance"                  ,
                      "rc_type_2":"Distance"                  ,
                      "adaptative":False                      , 
                      "save_format":".dcd"                    ,
                      "rmsGradient":0.1                       ,
                      "optimizer":"ConjugatedGradient"        ,
                      "MC_RC1":False                          ,
                      "MC_RC2":False                          ,
                      "log_frequency":50                      ,
                      "contour_lines":10                      ,
                      "nprocs":1                              ,
                      "show":False                             }
        





def parse_2D_scan_logfile (logfile):
    """ Function doc """
    
    data  =  open(logfile, 'r')
    lines =  data.readlines()
    xlist = []
    ylist = []
    zlist = []
    
    lastline = lines[-1].split()
    x_size = int(lastline[0])
    y_size = int(lastline[1])

    
    rows = y_size+1
    cols = x_size+1
     
    Z       = [[0]*cols for _ in range(rows)]
    RC1     = [[0]*cols for _ in range(rows)]
    RC2     = [[0]*cols for _ in range(rows)]
    
    #print (zlist)
    
    for line in lines[1:]:
        line2 = line.split()
        x = int(line2[0])
        y = int(line2[1])
        #print(x,y, line2[-1])
        
        
        
        Z[y][x]       = float(line2[-1]) 
        RC1[y][x]     = float(line2[-3]) 
        RC2[y][x]     = float(line2[-2]) 

    data = {
           'RC1': RC1,
           'RC2': RC2,
           'Z': Z
           }
    return data












    '''


    #zlist = [[None]*(x_size+1)]*(y_size+1)
    #
    #for line in lines[1:]:
    #    line2 = line.split()
    #    x = int(line2[0])
    #    y = int(line2[1])
    #    print(x,y, line2[-1])
    #    zlist[y][x] = float(line2[-1])        
    
    #print (zlist)
    
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors


    
    #self.Z = [[56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], 
    #          [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036], [56.34677450550225, 58.942831334607035, 73.21370758385456, 84.50838319960894, 84.47451804763841, 58.662963058639434, 38.687759070060565, 16.784363008395303, 33.923242461554764, 49.60649282127997, 45.99324563978007, 41.97248472983483, 40.403898982011015, 41.64380281602644, 44.31448492241907, 26.295264548542036]]
    self.Y  = range(0,len(self.Z))
    self.X  = range(0,len(self.Z[0]))
    
    
    #N = 10 
    #X, Y = np.mgrid[-3:3:complex(0, N), -2:2:complex(0, N)]
    #Z1 = np.exp(-X**2 - Y**2)
    #Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2)
    #Z = (Z1 - Z2) * 2
    #print(Z)
    
    #Y  = range(0,len(zlist))
    #X  = range(0,len(zlist[0]))
    #Z  = zlist 
    
    fig, ax = plt.subplots(1, 1)
    #pcm = ax[0].pcolormesh(X, Y, Z, norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=-1.0, vmax=1.0, ), cmap='RdBu_r', shading='nearest')
    #fig.colorbar(pcm, ax=ax[0], extend='both')

    pcm = ax.pcolormesh(self.X, self.Y, self.Z, cmap='jet', vmin=np.min(self.Z), shading='gouraud')
    fig.colorbar(pcm, ax=ax)#, extend='both')

    plt.show()
    '''

def main(args):
    
    win = PotentialEnergyAnalysisWindow()
    win.OpenWindow()
    Gtk.main()
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
