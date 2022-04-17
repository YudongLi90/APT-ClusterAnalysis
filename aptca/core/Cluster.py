"""
Container for atom cluster
"""

from platform import processor
from typing import List
import numpy as np
from scipy.spatial import Delaunay,ConvexHull
from sklearn import cluster
from aptca.core.APTData import APTData

from aptca.core.RRNG import RRNG
import multiprocessing
from aptca.utils.logging import log
from time import time
import datetime
from pathos.multiprocessing import ProcessingPool as Pool

class Cluster:

    ALL_POS = None
    ALL_MASS = None
    IONS_IN_CLUSTER = None
    IONS_TO_BE_ADDED = None
    POS_TO_BE_ADDED = None
    MASS_TO_BE_ADDED = None

    def __init__(self,
                label: int = None,
                pos :np.ndarray=None, 
                mass:np.ndarray=None,
                probabilities: np.ndarray=None,
                persistence:float = None,
                RRNG : RRNG = None
                ):
        self._label = label
        self.pos = pos
        self.mass = mass
        self.ions = None # ion name corresponding to mass
        self.probabilities = probabilities
        self._persistence = persistence
        self._RRNG = RRNG
        self.possible_ions = RRNG.ions # unique ions
        self.composition = None
        self.cluster_center_original = self._calc_cluster_center()
        self._get_cluster_polyhedron()
        self._calc_cluster_volume()
        self.atom_num_density_original = self._calc_atom_num_density()
    
    def _calc_cluster_center(self):
        centerX = np.mean(self.pos[:,0])
        centerY = np.mean(self.pos[:,1])
        centerZ = np.mean(self.pos[:,2])
        return np.array([centerX,centerY,centerZ])

    def _get_cluster_polyhedron(self):
        self._polyhedron = Delaunay(self.pos)
    
    def in_cluster(self,pos):
        """
        pos can be ndarray or single point(...,dim)
        """
        return self._polyhedron.find_simplex(pos) >= 0

    def _calc_cluster_volume(self):
        hull = ConvexHull(self.pos)
        self._volume = hull.volume
    
    def _calc_atom_num_density(self):
        """
        atom number density calculation, unit (count/nm3)
        """
        return self.pos.shape[0]/self._volume

    def get_cluster_simplified_info(self):
        """
        helper function to return cluster info for simplified data output
        """
        #label persistancy MinProb MeanProb MaxProb XClusCenter YClusCenter ZClusCenter ClusterSize 
        return [self.label,self.persistence,min(self.probabilities),np.mean(self.probabilities),max(self.probabilities),
                self.cluster_center_original[0],self.cluster_center_original[1],self.cluster_center_original[2],len(self.mass),
                ]

    def get_cluster_simplified_info_with_volume(self):
        #label persistancy MinProb MeanProb MaxProb XClusCenter YClusCenter ZClusCenter ClusterSize ClusterVolume ClusterAtom#Density
        return [self.label,self.persistence,min(self.probabilities),np.mean(self.probabilities),max(self.probabilities),
                self.cluster_center_original[0],self.cluster_center_original[1],self.cluster_center_original[2],len(self.mass),
                self._volume, self.atom_num_density_original
                ]

    
    def add_atom(self,
                pos: np.ndarray = None,
                mass: float = None,
                ):
        """
        add atom to existing cluster
        
        """
        if self.in_cluster(pos):
            self.pos = np.array([self.pos,pos])
            self.mass = np.array([self.mass,mass])          
    
    def add_atoms_bulk(self,pos_array, mass_list):
        index = self.in_cluster(pos_array)
        self.pos = np.concatenate((self.pos,pos_array[index]),axis=0)
        self.mass = np.concatenate((self.mass, mass_list[index]),axis=0)
        
        
    def add_missing_atoms(self):
        self.add_atoms_bulk(Cluster.POS_TO_BE_ADDED,Cluster.MASS_TO_BE_ADDED)

    def compute_composition(self):
        """
        composition is count of each type of atom (ion) in a cluster
        """
        composition = []
        for ion in Cluster.IONS_IN_CLUSTER:
            count = 0
            ranges = self._RRNG.get_range_for_ion(ion)
            for range_ion in ranges:
                index = np.where(np.logical_and(self.mass>=range_ion[0], self.mass<=range_ion[1]))[0]
                count += len(self.mass[index])
                #TODO update ion?
            #composition.append(count/len(self.mass))
            # Yudong: Apr.16 2022. Change composition to counts
            composition.append(count)
        self.composition = np.array(composition)
    
    def get_composition_from_ion_lst(self,ions : List[str]=None):
        composition = []
        if isinstance(ions,str) or isinstance(ions,list):
            # add missing ions to cluster
            if self.composition is None:
                self.compute_composition()
            elif len(self.composition) < len(Cluster.IONS_IN_CLUSTER):
                log.info("Computing compositions due to added ions")
                self.compute_composition()
            elif len(self.composition) > len(Cluster.IONS_IN_CLUSTER):
                raise ValueError("Please doubleckec. composition has more value than IONS_IN_CLUSTER")
            elif len(self.composition) == len(Cluster.IONS_IN_CLUSTER):
                pass
            else:
                log.warn("Something is wrong. Check get_composition_from_ion_lst fucntion")

            ions_in_cluster = np.array(Cluster.IONS_IN_CLUSTER)
            for ion in ions:
                #composition.append(self.composition[np.where(ions_in_cluster==ion)[0]])
                # Yudong 2022Apr16. should only have 1 value in composition. Not checking for now
                composition.append(self.composition[np.where(ions_in_cluster==ion)[0]][0])         
        else:
            log.warn("Please provice correct ions list (as a list of strings).")
        return dict(ion=ions,composition=composition)
    
    # def get_composition_from_ion_lst(self,ions : List[str]=None):
    #     composition = []
    #     if isinstance(ions,str):
    #         composition = 1
    #     elif isinstance(ions,list):
    #         if len(ions) == 1:
    #             composition = 1
    #         else:
    #             # add missing ions to cluster
    #             if self.composition is None:
    #                 self.compute_composition()
    #             elif len(self.composition) < len(Cluster.IONS_IN_CLUSTER):
    #                 log.info("Computing compositions due to added ions")
    #                 self.compute_composition()
    #             elif len(self.composition) > len(Cluster.IONS_IN_CLUSTER):
    #                 raise ValueError("Please doubleckec. composition has more value than IONS_IN_CLUSTER")
    #             elif len(self.composition) == len(Cluster.IONS_IN_CLUSTER):
    #                 pass
    #             else:
    #                 log.warn("Something is wrong. Check get_composition_from_ion_lst fucntion")

    #             ions_in_cluster = np.array(Cluster.IONS_IN_CLUSTER)
    #             for ion in ions:
    #                 composition.append(self.composition[np.where(ions_in_cluster==ion)[0]])
                    
    #     else:
    #         log.warn("Please provice correct ions list (as a list of strings).")
    #     return dict(ion=ions,composition=composition)

    @property
    def label(self):
        return self._label

    @property
    def persistence(self):
        return self._persistence

    @property
    def cluster_size(self):
        return len(self.mass)

    @classmethod
    def set_all_pos_mass(cls,aptdata: APTData):
        cls.ALL_POS = aptdata.pos
        cls.ALL_MASS = aptdata.mass
    
    @classmethod
    def update_ions_in_cluster(cls,ions:List[str]):
        if cls.IONS_IN_CLUSTER == None:
            if isinstance(ions,str):
                cls.IONS_IN_CLUSTER = ions.split()
            elif isinstance(ions,list):
                cls.IONS_IN_CLUSTER = ions
            else:
                raise ValueError("Incorrect type received for ions (should be list).")
        else:
            for ion in ions:
                if ion not in cls.IONS_IN_CLUSTER:
                    cls.IONS_IN_CLUSTER.append(ion)
    
    @classmethod
    def get_ions_to_add_in_cluster(cls,ions: List[str]):
        cls.IONS_TO_BE_ADDED = None # reset this
        ions_to_be_added = []
        for ion in ions:
            if ion not in cls.IONS_IN_CLUSTER:
                ions_to_be_added.append(ion)
        cls.IONS_TO_BE_ADDED = ions_to_be_added
    
    @classmethod
    def after_ions_added(cls):
        # call this after all clusters are updated.
        cls.update_ions_in_cluster(cls.IONS_TO_BE_ADDED)
        cls.IONS_TO_BE_ADDED = None
        cls.POS_TO_BE_ADDED = None
        cls.MASS_TO_BE_ADDED = None
        
    
    @classmethod
    def set_pos_mass_to_be_added(cls,pos,mass):
        cls.POS_TO_BE_ADDED = pos
        cls.MASS_TO_BE_ADDED = mass

        



    
    
                
    

            

