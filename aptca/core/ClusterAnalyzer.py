"""
Perform HDBscan cluster analysis
"""

from pathos.multiprocessing import ProcessingPool as Pool
from typing import List
import hdbscan as hdb
import numpy as np
import pandas as pd

from aptca.core.APTData import APTData
from aptca.core.Cluster import Cluster
from aptca.core.RRNG import RRNG
from aptca.utils.visualizer import visualization
from aptca.utils.logging import log

from time import time
import datetime


def worker_add_atoms_to_cluster(cluster:Cluster):
    cluster.add_missing_atoms()

class ClusterAnalyzer:

    def __init__(self,min_cluster_size,min_samples,cluster_selection_method='eom',n_jobs:int = 1,
                                                approx_min_span_tree=False):
        self.n_jobs = n_jobs
        self.labels = None
        self.pos = None
        self.ions_for_fitting = None
        self.all_cluster_list = None
        self.cluster_list = None
        self.all_cluster_simplified_info = None
        self.filtered_cluster_simplified_info = None
        self._clusterer = hdb.HDBSCAN(min_cluster_size=min_cluster_size,min_samples=min_samples,cluster_selection_method=cluster_selection_method,
                                                approx_min_span_tree=approx_min_span_tree)
    
    def fit(self,aptdata:APTData,ion:str or List(str) = None):
        self.aptdata = aptdata
        self.ions_for_fitting = ion
        if ion == None:
            log.warn("No Ion selected for cluster analysis. Fitting for cluster using all atoms. May take a long time")
            self._full_range_used = True
            self.pos = aptdata.pos
            self.mass = aptdata.mass
            self._clusterer.fit(self.pos)
        elif isinstance(ion,str):
            self._full_range_used = False
            log.info(f"Fitting for clusters using {ion} atom positions.")
            self.pos = aptdata.filter_ion(ion).pos
            self.mass = aptdata.filter_ion(ion).mass
            self._clusterer.fit(self.pos)
        elif isinstance(ion,list):
            self._full_range_used = False
            log.info(f"Fitting for clusters using {ion} atom positions.")
            self.pos = aptdata.filter_ions(ion).pos
            self.mass = aptdata.filter_ions(ion).mass
            self._clusterer.fit(self.pos)
        else:
            raise ValueError("Please provide corect ion format either a string or list")

        unique_labels = np.unique(self._clusterer.labels_)
        non_noise_unique_labels = unique_labels[unique_labels >= 0]
        self.all_cluster_list = []
        all_cluster_simplified_info = []
        #for label in unique_labels:
        # not considering noise from HDBScan result
        for label in non_noise_unique_labels:
            index = np.where(self._clusterer.labels_ == label)[0]
            pos = self.pos[index,:]
            mass = self.mass[index]
            probabilities = self._clusterer.probabilities_[index]
            persistence = self._clusterer.cluster_persistence_[label]
            cluster = Cluster(label = label,pos = pos,mass=mass,probabilities = probabilities,
                                persistence=persistence, RRNG = self.aptdata.RRNG)
            
            all_cluster_simplified_info.append(cluster.get_cluster_simplified_info())
            self.all_cluster_list.append(cluster)
        # Yudong 2022 Apr 17 : add cluster_lst to all_cluster_list to avoid error plot cluster center (less hassle)
        self.cluster_lst = self.all_cluster_list

        self.all_cluster_simplified_info = np.array(all_cluster_simplified_info)
        self._all_persistence = self._clusterer.cluster_persistence_
        self.labels = self._clusterer.labels_
        self.cluster_atom_index = np.where(self._clusterer.labels_ >= 0)
        self.noise_index = np.where(self._clusterer.labels_ == -1)

    def filter_clusters(self,persistenceThreshold,probabilityThreshold):
        unique_labels = list(set(self._clusterer.labels_))
        #self._qualifiedLabels = unique_labels[np.where(self._clusterer.cluster_persistence_>persistenceThreshold)[0]]
        self._qualifiedLabels = np.where(self._clusterer.cluster_persistence_>=persistenceThreshold)[0] # this returns the index is the same as labels
        self.cluster_list = []
        filtered_cluster_simplified_info = []
        for label in self._qualifiedLabels:
            index = np.where(np.logical_and(self._clusterer.labels_ == label,
                                                self._clusterer.probabilities_>=probabilityThreshold
                                                ))[0]
            pos = self.pos[index,:]
            mass = self.mass[index]
            probabilities = self._clusterer.probabilities_[index]
            persistence = self._clusterer.cluster_persistence_[label]
            cluster = Cluster(label = label,pos = pos,mass=mass,probabilities = probabilities,
                                persistence=persistence, RRNG = self.aptdata.RRNG)
            filtered_cluster_simplified_info.append(cluster.get_cluster_simplified_info_with_volume())
            self.cluster_list.append(cluster)
        
        self.filtered_cluster_simplified_info = np.array(filtered_cluster_simplified_info)

        self.cluster_atom_index = np.where(np.logical_and(np.isin(self._clusterer.labels_, self._qualifiedLabels),
                                                        self._clusterer.probabilities_>= probabilityThreshold))
        self.noise_index = np.where(np.logical_or(~np.isin(self._clusterer.labels_,self._qualifiedLabels),
                                                self._clusterer.probabilities_<probabilityThreshold
                                                ))[0]
        #self.noise_pos = self.pos[self.noise_index]
        #self.noise_mass = self.mass[self.noise_index]
        Cluster.set_all_pos_mass(self.aptdata)
        Cluster.update_ions_in_cluster(self.ions_for_fitting)

    def job_add_atoms_to_cluster(self,cluster:Cluster):
        cluster.add_atoms_bulk(self.aptdata.pos,self.aptdata.mass)
    
    def parallel_add_atoms_to_clusters(self,ions: List[str]=None):
        #prepare for cluster to add missing atoms (ions)
        Cluster.get_ions_to_add_in_cluster(ions)
        aptdata_to_add = self.aptdata.filter_ions(Cluster.IONS_TO_BE_ADDED)
        Cluster.set_pos_mass_to_be_added(aptdata_to_add.pos,aptdata_to_add.mass)
        # perform atom addtion to each cluster.
        cluster_list = self.cluster_list
        log.info("parallel")
        print("Start!")
        start = time()
        #pool = multiprocessing.Pool(processes = 6)
        pool = Pool(8)
        #pool.map(self.job_add_atoms_to_cluster,self.cluster_list)
        pool.map(worker_add_atoms_to_cluster,cluster_list)
        time_elapsed = time()-start
        log.info(f"Parallel took time : {str(datetime.timedelta(seconds=time_elapsed))}")
        Cluster.after_ions_added()


    def serial_add_atoms_to_clusters(self,ions: List[str]=None):
        #for cluster in self.cluster_list:
        #    pass
        #prepare for cluster to add missing atoms (ions)
        #log.info("Preparing info for adding atoms of interest to clusters.")
        Cluster.get_ions_to_add_in_cluster(ions)
        aptdata_to_add = self.aptdata.filter_ions(Cluster.IONS_TO_BE_ADDED)
        Cluster.set_pos_mass_to_be_added(aptdata_to_add.pos,aptdata_to_add.mass)

        log.info("Start cluster compositional analysis.")
        start = time()
        for i, cluster in enumerate(self.cluster_list):
            worker_add_atoms_to_cluster(cluster)
            if i%20 == 0:
                log.info(f"{100*i/len(self.cluster_list):0.2f}% finished.")
        time_elapsed = time()-start
        log.info(f"Time elapsed: {str(datetime.timedelta(seconds=time_elapsed))}")
        #print(f"Time elapsed: {str(datetime.timedelta(seconds=time_elapsed))}")
        Cluster.after_ions_added()

    
    def compute_cluster_composition(self,ions: List[str]=None):
        composition = []
        labels = []
        for cluster in self.cluster_list:
            labels.append(cluster.label)
            composition.append(cluster.get_composition_from_ion_lst(ions))
        
        composition_array = np.concatenate((np.array(labels).reshape(-1,1),np.array([d['composition'] for d in composition])),axis=1)
        column_lst =  composition[0]['ion']
        column_lst.insert(0,'label')
        df_composition = pd.DataFrame(data=composition_array,columns=column_lst)
        self.df_composition = df_composition
        return df_composition

    def save_composition_result(self,filename):
        """
        saves result from composition analysis.
        argument is the return from self.compute_cluster_composition
        """
        if self.df_composition is None:
            raise ValueError("Please provide data frame for compositional analysis or perform compositional analysis first.")
        else:
            log.info(f"Saving compositional analysis information to file {filename} ")
            try:
                self.df_composition.to_excel(filename)
            except Exception as e:
                log.warning("Error occured in saving.")
                log.warning(e)

    def save_all_cluster_simplified_info(self,filename_txt):
        if self.all_cluster_simplified_info is None:
            raise ValueError('Please do cluster analysis before saving cluster info.')
        else:
            log.info("Saving all cluster information from cluster analyzer.")
            try:
                np.savetxt(filename_txt,self.all_cluster_simplified_info,fmt=['%6i','%10.3f','%10.2f','%10.2f','%10.2f','%10.5f','%10.5f','%10.5f','%6i'],delimiter='\t',comments='')
            except Exception as e:
                log.warning("Error occured in saving.")
                log.warning(e)
    
    def save_filtered_cluster_simplified_info(self,filename_txt):
        if self.filtered_cluster_simplified_info is None:
            raise ValueError('Please do cluster filter analysis before saving cluster info.')
        else:
            log.info("Saving filtered cluster information from cluster analyzer.")
            try:
                np.savetxt(filename_txt,self.filtered_cluster_simplified_info,fmt=['%6i','%10.3f','%10.2f','%10.2f','%10.2f','%10.5f','%10.5f','%10.5f','%6i','%10.5f','%10.5f'],delimiter='\t',comments='')
            except Exception as e:
                log.warning("Error occured in saving.")
                log.warning(e)
               

    def show_all_atoms(self,background=True):
        visualization(pos=self.pos,label_lst=self._clusterer.labels_,background=background)

    @property
    def all_persistence(self):
        if hasattr(self,'_all_persistence'):
            return self._all_persistence
        else:
            raise ValueError("Fit for clusters first.")
        

if __name__ == '__main__':
    aptData = APTData()
    apt_file = "./APTdata/43334-1.6-7.apt"
    #pos_file = "./APTdata/43441-1.58-7.6.POS"
    rrng_file = "./APTdata/43334.RRNG"
    aptData.load_data_fromfile(apt_file=apt_file,rrng_file=rrng_file)
    clusterAN = ClusterAnalyzer(min_cluster_size=15,min_samples=3)
    clusterAN.fit(aptData,'Mg1')
    clusterAN.filter_clusters(persistenceThreshold=0.042,probabilityThreshold=0.8)
    #clusterAN.parallel_add_atoms_to_clusters(['Cu1','Si1'])
    clusterAN.serial_add_atoms_to_clusters(['Cu1','Si1','Al1'])
    #cluster = clusterAN.cluster_list[0]
    #cluster.paralle_add_atoms_bulk(aptData.pos,aptData.mass)  

        

    


