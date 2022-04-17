from itertools import cycle
from typing import List
from PyQt5.QtCore import QObject, pyqtSlot, pyqtSignal
import pyqtgraph.opengl as gl
from pyqtgraph import Vector

import numpy as np
from aptca.GUI.mplwidget import MplWidget
from aptca.core.ClusterAnalyzer import ClusterAnalyzer
from aptca.utils.logging import log
from aptca.utils import ColorTable

def generate_colors(pos:np.ndarray,labels:List) -> np.ndarray:
    color_table = ColorTable.gen_color_table(style='float')
    color_table = np.delete(color_table, [0], axis=0)
    cycol = cycle(color_table)
    cluster_color = dict()
    unique_cluster_id = np.unique(labels)

    for unique_id in unique_cluster_id:
        if unique_id == -1:
            cluster_color[unique_id] = [0.0, 0.0, 0.0]
        else:
            cluster_color[unique_id] = next(cycol)

    color = np.empty((pos.shape[0], 4))
    color[:, 3] = 0.1
    for idx, c_id in enumerate(labels):
        color[idx, 0:3] = cluster_color[c_id]
    return color


class AtomVisualizationWorker(QObject):
    data_ready = pyqtSignal(object)

    def __init__(self,pos,GLView,mplWidget: MplWidget=None,size=1,noise_color=[0.0,0.0,0.0,0.5],pxMode=True):
        super().__init__()
        self.view = GLView
        self.mplWidget = mplWidget
        self.fig = mplWidget.canvas.fig
        self.size = size
        self.noise_color=noise_color
        self.pxMode = pxMode
        #data = pos.astype(np.float32)
        #data[:, 2] = -data[:, 2] # reverse the z axis to show needle right.
        # this only shows atoms filtered for cluster analysis 
        #self.pos = data
        self.pos = pos
        self.pos_original = pos
        self.clusterAnalyzer = None
        self.cluster_pos = None
        self.cluster_center = None
        self.x_range, self.y_range, self.z_range, self.center = self.get_data_range_and_center()
        self.x_location = self.x_range[0]
        self.y_location = self.y_range[0]
        self.z_location = self.z_range[0]
        # default plot parameters
        self.atom_size = 1
        self.cluster_size = 3
        self.cluster_center_size = 5
        self.slice_thickness = -1
        self.location = -1
        self.by_axis = -1
        self.set_distance_once = True
    
    def set_clusterAnalyzer(self,clusterAnalyzer: ClusterAnalyzer):
        self.clusterAnalyzer = clusterAnalyzer 
        self.cluster_pos = self.pos[clusterAnalyzer.cluster_atom_index]

    def get_data_range_and_center(self) -> List:
        minX = min(self.pos[:,0])
        maxX = max(self.pos[:,0])
        minY = min(self.pos[:,1])
        maxY = max(self.pos[:,1])
        minZ = min(self.pos[:,2])
        maxZ = max(self.pos[:,2])
        centerX = minX + (minX+maxX)/2
        centerY = minY + (minY+maxY)/2
        centerZ = minZ + (minZ+maxZ)/2
        return [minX, maxX], [minY, maxY] , [minZ, maxZ], [centerX, centerY, centerZ]

    def slice_data(self,pos,slice_thickness,location,by_axis) -> np.ndarray:
        self.slice_thickness = slice_thickness
        thickness = float(self.slice_thickness)
        if thickness <= 0:
            log.info("Non positive thickness. Showing all atoms.")
            #return np.linspace(0,pos.shape[0]-1,pos.shape[0],dtype=int)
            return None
        else:
            range_to_show = [location - thickness/2, location + thickness/2]
            index = np.where(np.logical_and(pos[:,by_axis]>=range_to_show[0],pos[:,by_axis]<=range_to_show[1]))[0]
            return index
    
    def show_atoms(self,size=1, slice_index = None):
        # this only shows atoms filtered for cluster analysis 
        self.view.clear()
        if slice_index is None:
            sp1 = gl.GLScatterPlotItem(pos=self.pos, size=size, color=self.noise_color, pxMode=self.pxMode)
            if self.mplWidget is not None:
                    fig = self.mplWidget.canvas.fig
                    fig.clf()
                    ax = fig.add_subplot(projection='3d')
                    ax.scatter(self.pos[:,0],self.pos[:,1],self.pos[:,2],s=size,color=self.noise_color)
                    ax.set_xlabel("X (nm)")
                    ax.set_ylabel("Y (nm)")
                    ax.set_zlabel("Z (nm)")
                    self.fig.tight_layout()
                    self.mplWidget.canvas.draw()
        else:
            sp1 = gl.GLScatterPlotItem(pos=self.pos[slice_index,:], size=size, color=self.noise_color, pxMode=self.pxMode)
            if self.mplWidget is not None:
                fig = self.mplWidget.canvas.fig
                fig.clf()
                pos = self.pos[slice_index,:]
                ax = fig.add_subplot(projection='3d')
                ax.scatter(pos[:,0],pos[:,1],pos[:,2],s=size,color=self.noise_color)
                ax.set_xlabel("X (nm)")
                ax.set_ylabel("Y (nm)")
                ax.set_zlabel("Z (nm)")
                self.fig.tight_layout()
                self.mplWidget.canvas.draw()
        sp1.setGLOptions('opaque')
        self.view.setBackgroundColor('w')
        if self.set_distance_once:
            self.view.opts['distance'] = 3000
            self.view.opts['center'] = Vector(self.center)
            self.set_distance_once = False

        #self.view.opts['fov'] = 1
        self.view.addItem(sp1)

    def show_clusters(self, clusterAnalyzer: ClusterAnalyzer,cluster_size, noise_size,cluster_center_size,
                        cluster_slice_index=None,noise_slice_index=None,
                        cluster_center_slice_index = None):
        # this only shows atoms filtered for cluster analysis
        self.view.clear()
        self.clusterAnalyzer = clusterAnalyzer
        self.cluster_pos = self.pos[clusterAnalyzer.cluster_atom_index]
        self.noise_pos = self.pos[clusterAnalyzer.noise_index]
        self.cluster_center = np.array([cluster._calc_cluster_center().tolist() for cluster in clusterAnalyzer.cluster_lst])
        self.cluster_center_label = [cluster.label for cluster in clusterAnalyzer.cluster_lst]
        cluster_index = clusterAnalyzer.cluster_atom_index
        noise_index = clusterAnalyzer.noise_index
        cluster_pos = self.pos[cluster_index]
        cluster_labels = self.clusterAnalyzer.labels[cluster_index]
        noise_pos = self.pos[noise_index]
        cluster_color = generate_colors(cluster_pos,cluster_labels)
        cluster_center_color = generate_colors(self.cluster_center,self.cluster_center_label)

        if cluster_slice_index is None:
            #show all tip
            # show noise atoms
            if noise_size > 0:
                sp1 = gl.GLScatterPlotItem(pos=noise_pos,size=noise_size,color=self.noise_color,pxMode=self.pxMode)
                sp1.setGLOptions('translucent')
                self.view.addItem(sp1)
            if cluster_size > 0:
                sp2 = gl.GLScatterPlotItem(pos=cluster_pos,size=cluster_size,color=cluster_color,pxMode=self.pxMode)
                sp2.setGLOptions('opaque')
                self.view.addItem(sp2) 
            if cluster_center_size > 0:
                sp3 = gl.GLScatterPlotItem(pos=self.cluster_center,size=cluster_center_size,color=cluster_center_color,pxMode=self.pxMode)
                sp3.setGLOptions('opaque')
                self.view.addItem(sp3)    
        elif (cluster_slice_index is not None) and (noise_slice_index is None):
            pass
        else:
            if noise_size > 0:
                sp1 = gl.GLScatterPlotItem(pos=noise_pos[noise_slice_index],size=noise_size,color=self.noise_color,pxMode=self.pxMode)
                sp1.setGLOptions('translucent')
                self.view.addItem(sp1)
            if cluster_size > 0:    
                sp2 = gl.GLScatterPlotItem(pos=cluster_pos[cluster_slice_index],size=cluster_size,
                                        color=cluster_color[cluster_slice_index],pxMode=self.pxMode)
                sp2.setGLOptions('opaque')
                self.view.addItem(sp2) 
            if cluster_center_size > 0:                            
                sp3 = gl.GLScatterPlotItem(pos=self.cluster_center[cluster_center_slice_index],size=cluster_center_size,
                                        color=cluster_center_color[cluster_center_slice_index],pxMode=self.pxMode)          
                sp3.setGLOptions('opaque')
                self.view.addItem(sp3)
            
        self.view.setBackgroundColor('w')



    def get_range_by_axis(self,axis):
        if axis == 0: # x axis
            return self.x_range
        elif axis == 1:
            return self.y_range
        elif axis == 2:
            return self.z_range
        else:
            return [0,100]

    def update_plot_parameters(self,updateview= True, **kwargs):
        allowed_keys = {'atom_size', 'cluster_size','cluster_center_size', 'slice_thickness','location','by_axis','clusterAnalyzer'}
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in allowed_keys)
        if updateview:
            self.update_view()


    def update_view(self):
        self.view.clear()
        if self.clusterAnalyzer is None:
            #TODO
            #self.show_atoms()
            slice_index = self.slice_data(self.pos,self.slice_thickness,self.location,self.by_axis) 
            self.show_atoms(size = self.atom_size, slice_index=slice_index)
        else:
            cluster_slice_index = self.slice_data(self.cluster_pos,self.slice_thickness,self.location,self.by_axis)
            noise_slice_index = self.slice_data(self.noise_pos,self.slice_thickness,self.location,self.by_axis)
            cluster_center_slice_index = self.slice_data(self.cluster_center,self.slice_thickness,self.location,self.by_axis)
            self.show_clusters(clusterAnalyzer=self.clusterAnalyzer,cluster_size = self.cluster_size,
                                noise_size=self.atom_size,cluster_center_size=self.cluster_center_size, cluster_slice_index=cluster_slice_index,
                                noise_slice_index=noise_slice_index,cluster_center_slice_index = cluster_center_slice_index)

            




