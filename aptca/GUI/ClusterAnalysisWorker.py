
from typing import List
from PyQt5.QtCore import QObject, QThread, pyqtSlot, pyqtSignal
from aptca.core.ClusterAnalyzer import ClusterAnalyzer
from aptca.core.APTData import APTData
from aptca.utils.logging import log
import pandas as pd

class ClusterAnalysisWorker(QObject):
    cluster_ready = pyqtSignal(object)
    composition_ready = pyqtSignal(object)
    finished = pyqtSignal()
    def __init__(self,aptdata:APTData,ions_selected: List[str],minClusterSize:float,minSamples:float):
        super().__init__()
        log.info(f"Cluster Analyzer initialized with minimum cluster size: {minClusterSize}, min samples: {minSamples}.")
        self.clusterAN = ClusterAnalyzer(min_cluster_size=minClusterSize,min_samples=minSamples)
        self.aptdata = aptdata
        self.ions_selected = ions_selected
        self.ions_for_compositional_analysis = None
        self.readyforcompositionalAnalysis = False
    
    @pyqtSlot(name="fit_clusters")
    def fit_clusters(self):
        try:
            self.clusterAN.fit(self.aptdata,self.ions_selected)
            self.cluster_ready.emit(self.clusterAN)
            self.readyforcompositionalAnalysis = True
        except Exception as e:
            log.warn(e)
            self.cluster_ready.emit(None)
        finally:
            self.finished.emit()
    
    #@pyqtSlot(name="filter_clusters")
    def filter_clusters(self,persistenceThreshold, probabilityThreshold):
        self.clusterAN.filter_clusters(persistenceThreshold=persistenceThreshold,probabilityThreshold=probabilityThreshold)
    
    def set_ions_for_composition_analysis(self,ion_lst: List[str]):
        if isinstance(ion_lst,List):
            self.ions_for_compositional_analysis = ion_lst
            log.info(f"{ion_lst} selected for compositional analysis.")
        else:
            log.warning("Check the implementation for selecting ions. Wrong type of arguments passed to set_ions_for_composition_analysis.")
    
    @pyqtSlot(name="compute_composition")
    def compute_composition(self):
        if self.readyforcompositionalAnalysis:
            if self.ions_for_compositional_analysis is None:
                log.warning("Please select ions for compositional analysis first!")
                self.cluster_ready.emit(None)
            else:
                log.info(f"Starting to perform compositional analysis for ions: {self.ions_for_compositional_analysis}.")
                self.clusterAN.serial_add_atoms_to_clusters(self.ions_for_compositional_analysis)
                composition_dict = self.clusterAN.compute_cluster_composition(self.ions_for_compositional_analysis)
                #df = pd.DataFrame(composition_dict)
                #df.to_pickle() get path from main.
                self.composition_ready.emit(self.clusterAN)
        else:
            log.warning("Please first perform cluster analysis!")
            self.cluster_ready.emit(None)
        self.finished.emit()




    
    





