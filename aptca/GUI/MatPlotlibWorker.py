from xml.sax.handler import all_properties
import matplotlib.pyplot as plt
from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot
import numpy as np
from sklearn import cluster
from aptca.core.APTData import APTData
from aptca.core.ClusterAnalyzer import ClusterAnalyzer
from aptca.GUI.mplwidget import MplWidget

class MatPlotlibWorker:
    """
    just a normal class since you cannot draw figure in embeded GUI in a seperate thread
    """
    def __init__(self,aptdata: APTData, mplWidget: MplWidget):
        super().__init__()
        self.aptdata = aptdata
        self.mplWidget = mplWidget
        self.fig = mplWidget.canvas.fig
    
    @pyqtSlot(name="plot_atomDensity")
    def plot_atomDensity(self):
        self.fig.clf() # clear before plot
        distance_axial = self.aptdata.atomcloud.distance_axial
        density_axial = self.aptdata.atomcloud.avg_density_axial
        distance_radial = self.aptdata.atomcloud.distance_radial
        density_radial = self.aptdata.atomcloud.avg_density_radial
        Density = self.aptdata.atomcloud.Density
        axial_range = self.aptdata.atomcloud.axial_range
        radial_range = self.aptdata.atomcloud.radial_range
        axial_range_1_8 = self.aptdata.atomcloud.axial_density_Range_1_8
        radial_range_1_8 = self.aptdata.atomcloud.radial_density_Range_1_8
        ax0 = self.fig.add_subplot(3,1,1)
        #fig.suptitle("Atom Density Distributions",fontsize=20)
        # plot axial density distribution
        ax0.plot(distance_axial,density_axial,color='k',linewidth=2)
        ax0.set_ylabel('Density ($\mathregular{at/nm^3}$)')
        ax0.set_xlabel('Depth (nm)')
        ax0.set_title(f'Axial Density Distribution, range: {axial_range:0.4f}',fontsize=16)
        # radial
        ax1 = self.fig.add_subplot(3,1,2)
        ax1.plot(distance_radial,density_radial,color='k',linewidth=2)
        ax1.set_ylabel('Density ($\mathregular{at/nm^3}$)')
        ax1.set_xlabel('Radius (nm)')
        ax1.set_title(f'Radial Density Distribution,range:{radial_range:0.4f}',fontsize=16)
        # histogram
        ax3 = self.fig.add_subplot(3,1,3)
        ax3.hist(Density,bins=100)
        ax3.set_ylabel('Counts')
        ax3.set_xlabel('Density ($\mathregular{at/nm^3}$)')
        ax3.set_title(f'Density Histogram, old ranges, axial: {axial_range_1_8:0.4f}, radial: {radial_range_1_8:0.4f}',fontsize=16)
        self.fig.tight_layout()
        self.mplWidget.canvas.draw()
        #self.finished.emit()

    def plot_clustersize_v_persistence(self,clusterAnalyzer: ClusterAnalyzer):
        self.fig.clf() # clear before plot
        clustersize = [cluster.cluster_size for cluster in clusterAnalyzer.all_cluster_list]
        persistence = [cluster.persistence for cluster in clusterAnalyzer.all_cluster_list]
        ax = self.fig.add_subplot(1,1,1)
        ax.scatter(persistence, clustersize,color='k')
        ax.set_ylabel("Cluster size (# atoms)")
        ax.set_xlabel("Persistence")
        self.fig.tight_layout()
        self.mplWidget.canvas.draw()

    def plot_cumulativeNumberClusters_v_persistence(self,clusterAnalyzer:ClusterAnalyzer):
        self.fig.clf() # clear before plot
        persistence = np.array(clusterAnalyzer.all_persistence)
        persistence.sort()
        descending_p = persistence[::-1]
        cluster_counts = np.ones(len(descending_p))
        cumulative = np.cumsum(cluster_counts)
        
        ax = self.fig.add_subplot(1,1,1)
        ax.scatter(descending_p, cumulative,s=0.5,color='k')
        ax.set_ylabel("Cluster size (# atoms)")
        ax.set_xlabel("Persistence")
        self.fig.tight_layout()
        self.mplWidget.canvas.draw()

    def plot_atom_map_3D(self):
        raise NotImplementedError()



