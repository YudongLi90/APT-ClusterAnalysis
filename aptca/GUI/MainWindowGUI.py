import os
from PyQt5.QtWidgets import QApplication, QWidget, QMainWindow, QFileDialog, QGridLayout
from PyQt5 import uic
from PyQt5.QtCore import pyqtSignal, pyqtSlot, Qt, QThread, QObject, QPointF
import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from aptca.core.ClusterAnalyzer import ClusterAnalyzer

import ui.MainWindow as MainWindow

import sys
import requests

from aptca.core.APTData import APTData
from aptca.utils.visualizer import visualize_atom_cloud
from aptca.utils.logging import log
from aptca.GUI.AtomVisualizationWorker import AtomVisualizationWorker
from aptca.GUI.QtLogger import QtLogger
from aptca.GUI.GUIhelper import update_slider_and_label_range
from aptca.GUI.FileIOWorker import DataLoader
from aptca.GUI.ClusterAnalysisWorker import ClusterAnalysisWorker
from aptca.GUI.AtomDensityWorker import AtomDensityWorker
from aptca.GUI.MatPlotlibWorker import MatPlotlibWorker

class MainWindow(QMainWindow, MainWindow.Ui_MainWindow):

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.setupUi(self)
        self.aptData = None
        self.data_fname = None
        self.ions_selected_clusterA = None
        self.ions_selected_compositionalAnalysis = None
        self.selectedAtomCloud = None
        self.sliceThickness = -1
        self.sliceLocation = -1
        self.slicePerpendicularAxis = None # 0 for x-axis, 1 for y-axis, 2 for z-axis
        self.range = None
        self.atom_visualizer = None
        self.clusterAnalyzerWorker = None
        self.mplWorker = None
        self.minClusterSize = 15
        self.minSamples = 3
        self.persistenceThreshold = 0.042
        self.probabilityThreshold = 0.8
        self.atom_size = 1
        self.cluster_size = 3
        self.cluster_center_size = 5
        self.path = os.getcwd()
        self.log = QtLogger(self.plainTextEdit_log).get_logger()
        self.pushButton_openAPT.clicked.connect(self.button_openAPT)
        self.pushButton_openRRNG.clicked.connect(self.button_openRRNG)
        self.pushButton_loadData.clicked.connect(self.button_loadData)
        self.pushButton_showAtomDensity.clicked.connect(self.compute_atomDensity)
        self.listWidget_ionSelection.itemClicked.connect(self.show_selection)
        self.listWidget_ionSelectionForCompositonalAnalysis.itemClicked.connect(self.show_ion_selected_for_compositionAnalysis)
        self.pushButton_filterIonClusterA.clicked.connect(self.button_filterIon)
        self.pushButton_showAtoms.clicked.connect(self.button_showAtoms)
        self.spinBox_sliceThickness.valueChanged.connect(self.set_sliceThickness)
        self.spinBox_sliceLocation.valueChanged.connect(self.set_sliceLocation)
        self.radioButton_xAxis.toggled.connect(self.set_sliceAxis_to_x)
        self.radioButton_yAxis.toggled.connect(self.set_sliceAxis_to_y)
        self.radioButton_zAxis.toggled.connect(self.set_sliceAxis_to_z)
        self.pushButton_showAtomDensity.clicked.connect(self.compute_atomDensity)
        # clusterAnalysis section
        self.spinBox_minClusterSize.valueChanged.connect(self.set_minClusterSize)
        self.spinBox_minSamples.valueChanged.connect(self.set_minSamples)
        self.lineEdit_persistenceThreshold.editingFinished.connect(self.set_persistenceThreshold)
        self.lineEdit_probabilityThreshold.editingFinished.connect(self.set_probabilityThreshold)
        self.pushButton_fitClusters.clicked.connect(self.fit_clusters)
        self.pushButton_filterClusters.clicked.connect(lambda: self.filter_clusters(self.persistenceThreshold,self.probabilityThreshold))
        #compositonal analysis
        self.pushButton_CompositionAnalysis.clicked.connect(self.compositional_analysis)
        # visualization
        self.doubleSpinBox_atomSize.valueChanged.connect(self.set_atomSize)
        self.doubleSpinBox_clusterAtomSize.valueChanged.connect(self.set_clusterSize)
        self.doubleSpinBox_clusterCentroidSize.valueChanged.connect(self.set_clusterCenterSize)
        # save files
        self.pushButton_saveAllClusterInfo.clicked.connect(self.save_allClusterInfo)
        self.pushButton_saveFilteredClustersInfo.clicked.connect(self.save_filteredClusterInfo)
        self.pushButton_saveCompositionalInfo.clicked.connect(self.save_compositionalInfo)
        # converter
        self.actionAPT_to_POS.triggered.connect(self.convert_apt_to_pos)

        #update dummy values for range
        update_slider_and_label_range(self.horizontalSlider_sliceLocation,
                                        self.label_minLocation,self.label_maxLocation,[0,100])
        self.horizontalSlider_sliceLocation.valueChanged.connect(self.update_sliceLocation)
    
    def update_path(self,fname:str):
        self.path = os.path.dirname(os.path.abspath(fname))

    @pyqtSlot(name="actionAPT_to_POS")
    def convert_apt_to_pos(self):
        if self.aptData is None:
            raise ValueError("please load apt file first!")
        else:
            fname = QFileDialog.getSaveFileName(self,'Save apt info to pos',os.path.join(self.path,self.data_fname+'.POS'), "Text (*.POS)")
            save_fname = fname[0]
            if save_fname == "":
                log.info("Not saving")
            else:
                log.info(f"Saving to file: {save_fname}")
                try:
                    pos = self.aptData.atomcloud.pos
                    mass = self.aptData.atomcloud.mass
                    data = np.concatenate((pos,mass.reshape(-1,1)),axis=1)
                    self.aptData.atomcloud.write_pos(save_fname,data)
                    self.update_path(save_fname)
                except Exception as e:
                    log.warning("Save all cluster info failed.")
                    log.warning(e)

    @pyqtSlot(name="pushButton_saveAllClusterInfo")
    def save_allClusterInfo(self):
        fname = QFileDialog.getSaveFileName(self,'Save all cluster info to txt',os.path.join(self.path,'AllClusterInfo.txt'), "Text (*.txt)")
        save_fname = fname[0]
        if save_fname == "":
            log.info("Not saving")
        else:
            log.info(f"Saving to file: {save_fname}")
            if self.clusterAnalyzer is None:
                raise ValueError("Please perform cluster analysis before saving")
            else:
                try:
                    self.clusterAnalyzer.save_all_cluster_simplified_info(save_fname)
                    self.update_path(save_fname)
                except Exception as e:
                    log.warning("Save all cluster info failed.")
                    log.warning(e)
    
    @pyqtSlot(name="pushButton_saveFilteredClusterInfo")
    def save_filteredClusterInfo(self):
        fname = QFileDialog.getSaveFileName(self,'Save filtered cluster info to txt',os.path.join(self.path,'FilteredClusterInfo.txt'), "Text (*.txt)")
        save_fname = fname[0]
        if save_fname == "":
            log.info("Not saving")
        else:
            log.info(f"Saving to file: {save_fname}")
            if self.clusterAnalyzer is None:
                raise ValueError("Please perform cluster analysis before saving")
            else:
                try:
                    self.clusterAnalyzer.save_filtered_cluster_simplified_info(save_fname)
                    self.update_path(save_fname)
                except Exception as e:
                    log.warning("Save filtered cluster info failed.")
                    log.warning(e)
    
    @pyqtSlot(name="pushButton_saveCompositionalInfo")
    def save_compositionalInfo(self):
        fname = QFileDialog.getSaveFileName(self,'Save compositional info to txt',os.path.join(self.path,'clusterCompositionalInfo.txt'), "Text (*.txt)")
        save_fname = fname[0]
        if save_fname == "":
            log.info("Not saving")
        else:
            log.info(f"Saving to file: {save_fname}")
            if self.clusterAnalyzer is None:
                raise ValueError("Please perform cluster analysis before saving")
            else:
                try:
                    self.clusterAnalyzer.save_composition_result(save_fname)
                    self.update_path(save_fname)
                except Exception as e:
                    log.warning("Save cluster compositional info failed.")
                    log.warning(e)

        

    @pyqtSlot(name="pushButton_openAPT")
    def button_openAPT(self):
        fname = QFileDialog.getOpenFileName(self,'Open APT/POS file', self.path,"Files (*.APT or *.POS)")
        self.lineEdit_pos_filename.setText(fname[0])
        self.update_path(fname[0])

    @pyqtSlot(name="pushButton_openRRNG")
    def button_openRRNG(self):
        fname = QFileDialog.getOpenFileName(self,'Open RRNG file', self.path,"Files (*.RRNG)")
        self.lineEdit_RRNG_filename.setText(fname[0])

    @pyqtSlot(name="pushButton_loadData")
    def button_loadData(self):
        apt_fname = self.lineEdit_pos_filename.text()
        self.data_fname = '.'.join(os.path.basename(apt_fname).split('.')[0:-1])
        RRNG_fname = self.lineEdit_RRNG_filename.text()
        
        self.log.info(f"Position file: {os.path.basename(apt_fname)}")
        self.log.info(f"Range file: {os.path.basename(RRNG_fname)}")
        self.log.info("loading data from files.")

        self.dataworker = DataLoader(apt_fname=apt_fname,RRNG_fname=RRNG_fname)
        self.thread_dataworker = QThread()
        self.dataworker.data_ready.connect(self.loadAPTData)
        self.dataworker.finished.connect(self.thread_dataworker.quit)
        self.dataworker.moveToThread(self.thread_dataworker)
        self.thread_dataworker.started.connect(self.dataworker.load_data)
        self.thread_dataworker.start()
        #self.aptData = APTData()
        #self.aptData.load_data_fromfile(apt_file=apt_fname,rrng_file=RRNG_fname)
        #for ion in self.aptData.ions:
        #    self.listWidget_ionSelection.addItem(ion)

    @pyqtSlot(object, name='loadAPTData')
    def loadAPTData(self,data:APTData):
        if data is None:
            self.log.warn("Data load failed.")
        else:
            self.aptData = data
            for ion in self.aptData.ions:
                self.listWidget_ionSelection.addItem(ion)
                self.listWidget_ionSelectionForCompositonalAnalysis.addItem(ion)
            
    
    @pyqtSlot(name="compute_atomDensity")
    def compute_atomDensity(self):
        if self.aptData is None:
            log.warning("Please load APT data and RNNG files first!")
        else:
            log.info("Starting atom density worker")
            self.atomdensityworker = AtomDensityWorker(self.aptData)
            self.thread_ADworker = QThread(parent=self)
            self.atomdensityworker.data_ready.connect(self.get_atomDensity)
            self.atomdensityworker.finished.connect(self.thread_dataworker.quit)
            self.atomdensityworker.moveToThread(self.thread_ADworker)
            self.thread_ADworker.started.connect(self.atomdensityworker.calc_atomdensity)
            self.thread_ADworker.start()
    
    @pyqtSlot(object,name="get_atomDensity")
    def get_atomDensity(self,data):
        if data is None:
            log.warn("Atom density calculation failed.")
        else:
            self.aptData = data # update the aptData after AtomDensity calculation.
            log.info("Atom density cacluation finished. Plot results.")
            self.plot_atomdensity()

    def plot_atomdensity(self):
        self.mplWorker = MatPlotlibWorker(self.aptData,self.gV_MplWidget)
        self.mplWorker.plot_atomDensity()


    def show_selection(self):
        items = self.listWidget_ionSelection.selectedItems()
        x = []
        for i in range(len(items)):
            x.append(str(self.listWidget_ionSelection.selectedItems()[i].text()))
        self.ions_selected_clusterA = x
        self.log.info(f"Ions selected for cluster analysis: {x}")
        
    def show_ion_selected_for_compositionAnalysis(self):
        items = self.listWidget_ionSelectionForCompositonalAnalysis.selectedItems()
        x = []
        for i in range(len(items)):
            x.append(str(self.listWidget_ionSelectionForCompositonalAnalysis.selectedItems()[i].text()))
        self.ions_selected_compositionalAnalysis = x
        self.log.info(f"Ions selected for composition analysis: {x}")

    #########
    # cluster analysis
    #######
    @pyqtSlot(name='set_minClusterSize')
    def set_minClusterSize(self):
        self.minClusterSize = self.spinBox_minClusterSize.value()
        log.info(f"User set minimum cluster size to {self.minClusterSize} atoms/cluster")
    
    @pyqtSlot(name='set_minSamples')
    def set_minSamples(self):
        self.minSamples = self.spinBox_minSamples.value()
        log.info(f"User set minimum samples to {self.minSamples}, type: {type(self.minSamples)}")

    @pyqtSlot(name='set_persistenceThreshold')
    def set_persistenceThreshold(self):
        self.persistenceThreshold = float(self.lineEdit_persistenceThreshold.text())
        log.info(f"User set persistence threshold to {self.persistenceThreshold}")

    @pyqtSlot(name='set_probabilityThreshold')
    def set_probabilityThreshold(self):
        self.probabilityThreshold = float(self.lineEdit_probabilityThreshold.text())
        log.info(f"User set probability threshold to {self.probabilityThreshold}")
 
    @pyqtSlot(name="pushButton_filterIon")
    def button_filterIon(self):
        if self.aptData != None:
            self.selectedAtomCloud = self.aptData.filter_ions(ions=self.ions_selected_clusterA)
        else:
            self.log.warn("Please load data file first!")
    
    @pyqtSlot(name="fit_clusters")
    def fit_clusters(self):
        if self.minClusterSize is None and self.minSamples is None:
            log.info("Please set cluster fitting parameters first: minClusterSize & minSamples")
        else:
            self.clusterAnalyzerWorker = ClusterAnalysisWorker(aptdata=self.aptData,ions_selected=self.ions_selected_clusterA,
                                                                minClusterSize=self.minClusterSize,minSamples=self.minSamples)
            self.thread_CAworker = QThread()
            self.clusterAnalyzerWorker.cluster_ready.connect(self.get_clusters)
            self.clusterAnalyzerWorker.finished.connect(self.thread_CAworker.quit)
            self.clusterAnalyzerWorker.moveToThread(self.thread_CAworker)
            self.thread_CAworker.started.connect(self.clusterAnalyzerWorker.fit_clusters)
            self.thread_CAworker.start()

    @pyqtSlot(object,name="get_clusters")
    def get_clusters(self,data: ClusterAnalyzer):
        if data is None:
            self.log.warn("Cluster Analysis failed.")
        else:
            self.clusterAnalyzer = data
            self.log.info("Cluster Analysis done. Examine results.")
            if self.mplWorker is None:
                self.mplWorker = MatPlotlibWorker(self.aptData,self.gV_MplWidget)
            self.mplWorker.plot_cumulativeNumberClusters_v_persistence(self.clusterAnalyzer)
            if self.atom_visualizer is not None:
                self.atom_visualizer.show_clusters(clusterAnalyzer = self.clusterAnalyzer,
                                                    cluster_size= self.cluster_size,
                                                    noise_size= self.atom_size,
                                                    cluster_center_size= self.cluster_center_size,
                                                    cluster_slice_index=None
                                                    )

    
    @pyqtSlot(name="filter_clusters")
    def filter_clusters(self,persistenceThreshold, probabilityThreshold):
        if self.clusterAnalyzerWorker is None:
            log.warning("Fit cluster first before filtering.")
        else:
            self.clusterAnalyzerWorker.filter_clusters(persistenceThreshold,probabilityThreshold)
            if self.atom_visualizer is not None:
                self.atom_visualizer.show_clusters(clusterAnalyzer = self.clusterAnalyzerWorker.clusterAN,
                                                    cluster_size= self.cluster_size,
                                                    noise_size= self.atom_size,
                                                    cluster_center_size= self.cluster_center_size,
                                                    cluster_slice_index=None
                                                    )

    @pyqtSlot(name="compositional_analysis")
    def compositional_analysis(self):
        if self.ions_selected_compositionalAnalysis is None:
            log.warning("Please select ions for compositional analysis first.")
        elif self.clusterAnalyzerWorker is None:
            log.warning("Please perform cluster analysis before compositional analysis.")
        else:
            self.clusterAnalyzerWorker.set_ions_for_composition_analysis(self.ions_selected_compositionalAnalysis)

            self.thread_CompAworker = QThread()
            self.clusterAnalyzerWorker.cluster_ready.connect(self.get_compositionAnalysis)
            self.clusterAnalyzerWorker.finished.connect(self.thread_CompAworker.quit)
            self.clusterAnalyzerWorker.moveToThread(self.thread_CompAworker)
            self.thread_CompAworker.started.connect(self.clusterAnalyzerWorker.compute_composition)
            self.thread_CompAworker.start()
    
    @pyqtSlot(object,name="get_compositionAnalysis")
    def get_compositionAnalysis(self,data: ClusterAnalyzer):
        if data is None:
            self.log.warn("Compositional Analysis failed.")
        else:
            self.clusterAnalyzer = data
            self.log.info("Compositonal Analysis done. Examine results.")
            
        
    # visualization section
    @pyqtSlot(name="set_atomSize")
    def set_atomSize(self):
        self.atom_size = self.doubleSpinBox_atomSize.value()
        log.info(f"User set minimum atom size to {self.atom_size} ")
        if self.atom_visualizer is not None:
            self.atom_visualizer.update_plot_parameters(atom_size = self.atom_size)
    
    @pyqtSlot(name="set_clusterSize")
    def set_clusterSize(self):
        self.cluster_size = self.doubleSpinBox_clusterAtomSize.value()
        log.info(f"User set minimum cluster atom size to {self.cluster_size} ")
        if self.atom_visualizer is not None:
            self.atom_visualizer.update_plot_parameters(cluster_size = self.cluster_size)
    
    @pyqtSlot(name="set_clusterCenterSize")
    def set_clusterCenterSize(self):
        self.cluster_center_size = self.doubleSpinBox_clusterCentroidSize.value()
        log.info(f"User set minimum cluster atom size to {self.cluster_size} ")
        if self.atom_visualizer is not None:
            self.atom_visualizer.update_plot_parameters(cluster_center_size = self.cluster_center_size)


    
    @pyqtSlot(name="set_sliceThickness")
    def set_sliceThickness(self):
        self.sliceThickness = float(self.spinBox_sliceThickness.value())
        if self.sliceThickness > 0:
            log.info(f"User set slice thickness to {self.sliceThickness} nm.")
            if self.atom_visualizer is not None:
                self.atom_visualizer.update_plot_parameters(slice_thickness = self.sliceThickness)
            ##TODO
            ## update graph to reflect change
        else:
            log.info(f"User choose to show all atoms with thichness < 0nm.")
            self.slicePerpendicularAxis = -1
            if self.atom_visualizer is not None:
                self.atom_visualizer.update_plot_parameters(slice_thickness = self.sliceThickness,by_axis = self.slicePerpendicularAxis)

    @pyqtSlot(name="set_sliceLocation")
    def set_sliceLocation(self):
        self.sliceLocation = float(self.spinBox_sliceLocation.value())
        log.info(f"User set slice location to {self.sliceLocation} nm .")
        if self.atom_visualizer is not None:
            if self.sliceThickness <= 0:
                log.info(f"However, I cannot slice the result due to: User choose to show all atoms with thichness < 0nm.")
                self.atom_visualizer.update_plot_parameters(updateview=False,location=self.sliceLocation)
            else:
                self.atom_visualizer.update_plot_parameters(location=self.sliceLocation)
    
    @pyqtSlot(name="set_sliceAxis_to_x")
    def set_sliceAxis_to_x(self):
        if self.radioButton_xAxis.isChecked():
            self.slicePerpendicularAxis = 0
            if self.range != None:
                range_x = self.range[0]
                if self.atom_visualizer is not None:
                    self.atom_visualizer.update_plot_parameters(by_axis = self.slicePerpendicularAxis)
            else:
                self.log.warn("Data not loaded. Play with these dummy values [0,100] for range.")
                range_x = [0,100]
            update_slider_and_label_range(self.horizontalSlider_sliceLocation,
                                        self.label_minLocation,self.label_maxLocation,range_x)
            log.info(f"User choose to slice data perpendicular to X-axis.")
    
    @pyqtSlot(name="set_sliceAxis_to_y")
    def set_sliceAxis_to_y(self):
        if self.radioButton_yAxis.isChecked():
            self.slicePerpendicularAxis = 1
            if self.range != None:
                range_y = self.range[1]
                if self.atom_visualizer is not None:
                    self.atom_visualizer.update_plot_parameters(by_axis = self.slicePerpendicularAxis)
            else:
                self.log.warn("Data not loaded. Play with these dummy values [0,100] for range.")
                range_y = [0,100]
            update_slider_and_label_range(self.horizontalSlider_sliceLocation,
                                        self.label_minLocation,self.label_maxLocation,range_y)
            log.info(f"User choose to slice data perpendicular to Y-axis.")

    @pyqtSlot(name="set_sliceAxis_to_z")
    def set_sliceAxis_to_z(self):
        if self.radioButton_zAxis.isChecked():
            self.slicePerpendicularAxis = 2
            if self.range != None:
                range_z = self.range[2]
                if self.atom_visualizer is not None:
                    self.atom_visualizer.update_plot_parameters(by_axis = self.slicePerpendicularAxis)
            else:
                self.log.warn("Data not loaded. Play with these dummy values [0,100] for range.")
                range_z = [0,100]
            update_slider_and_label_range(self.horizontalSlider_sliceLocation,
                                        self.label_minLocation,self.label_maxLocation,range_z)
            log.info(f"User choose to slice data perpendicular to Z-axis.")
       
    
    @pyqtSlot(name="update_sliceLocation")
    def update_sliceLocation(self):
        self.slicelocation = self.horizontalSlider_sliceLocation.value()
        self.spinBox_sliceLocation.setValue(self.slicelocation)
        if self.atom_visualizer != None:
            #self.atom_visualizer.update_slice()
            if self.atom_visualizer is not None:
                    self.atom_visualizer.update_plot_parameters(location = self.sliceLocation)
        #log.info(f"slicer location {self.slicelocation}.")


    @pyqtSlot(name="pushButton_showAtoms")
    def button_showAtoms(self):
        if self.selectedAtomCloud != None:
            #visualize
            #visualize_atom_cloud(self.selectedAtomCloud.pos,self.graphicsView,size=1)
            self.atom_visualizer = AtomVisualizationWorker(pos=self.selectedAtomCloud.pos,GLView = self.graphicsView,
                                                            mplWidget=self.gV_MplWidget,size=1)
            self.atom_visualizer.show_atoms()
            self.spinBox_sliceThickness.setValue(-1)
            self.spinBox_sliceLocation.setValue(-1)
            x_range = self.atom_visualizer.get_range_by_axis(0)
            y_range = self.atom_visualizer.get_range_by_axis(1)
            z_range = self.atom_visualizer.get_range_by_axis(2)
            self.range=[x_range,y_range,z_range]
            range_default = self.atom_visualizer.get_range_by_axis(self.slicePerpendicularAxis)
            update_slider_and_label_range(self.horizontalSlider_sliceLocation,
                                        self.label_minLocation,self.label_maxLocation,range_default)
        else:
            axis = gl.GLAxisItem()
            self.graphicsView.addItem(axis)
    



if __name__ == "__main__":

    app = QApplication(sys.argv)
    form = MainWindow()
    form.show()
    app.exec_()