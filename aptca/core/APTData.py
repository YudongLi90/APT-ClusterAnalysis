"""
functions to read/write APT data/result
"""

from multiprocessing.sharedctypes import Value
import sys
import re
import os
from typing import Any, List, Set, Tuple, Union
import numpy as np
import pandas as pd
from sys import getsizeof
from functools import lru_cache
import os
import struct
from enum import Enum
from warnings import warn
from aptca.core.AtomCloud import AtomCloud
from aptca.utils import validate
from aptca.utils.logging import log
from aptca.core.RRNG import RRNG


class APTData():

    def __init__(self,atomcloud:AtomCloud = None,RRNG: RRNG = None,):
        """
        filename can be absolute or relative path + name (including file extension)
        pickel_path is the location where you store pickeled data
        cache is useful if you repeatly analyze one experiment (saves on data reading time) (default True)
        """
        self._atomcloud = atomcloud
        self._RRNG = RRNG
        self._files = {'Data filename':None, 'RRNG filename':None}

    def load_data_fromfile(self,apt_file:str, rrng_file:str):
        self._atomcloud = AtomCloud()
        self._atomcloud.load_atomcloud_fromfile(apt_file)
        self._RRNG = RRNG()
        self._RRNG.load_rrng_fromfile(rrng_file)
        #self._atomcloud.map_mass_to_ion(self._RRNG) # this is slow. Not using this.

    @property
    def pos(self) -> np.ndarray:
        return self._atomcloud.pos
    
    @property
    def mass(self) -> np.ndarray:
        return self._atomcloud.mass
    
    @property
    def mass_random(self) -> np.ndarray:
        return self._atomcloud.mass_random

    @property
    def elements(self) -> np.ndarray:
        return self._RRNG.elements
    
    @property
    def rrng(self) -> np.ndarray:
        return self._RRNG.rrng
    
    @property
    def ions(self) -> np.ndarray:
        return self._RRNG.ions

    @property
    def atomcloud(self) -> AtomCloud:
        return self._atomcloud
    
    @property
    def RRNG(self) -> RRNG:
        return self._RRNG
    
    @classmethod
    def _pickel(cls,df:pd.DataFrame,pickelpath,filename):
        if not os.path.exists(pickelpath) and pickelpath:
            print(f"Pickel path does not exist, creating folder: {pickelpath}")
            os.mkdir(pickelpath)
        df.to_pickle(pickelpath+'/'+filename)
        print(f"Saved data to a pickel file for fast loading next time: {pickelpath+'/'+filename}")

    def find_ion_indexes(self,ion:str) -> np.ndarray:
        ranges = self._RRNG.get_range_for_ion(ion)
        index = []
        for i,range in enumerate(ranges):
            index.append(np.where(np.logical_and(self._atomcloud._mass>=range[0], self._atomcloud._mass<=range[1]))[0])
        indexes = np.concatenate(index,axis=0)
        return indexes
    
    def filter_ion(self,ion:str) -> AtomCloud:
        """
        filter APT data based on ion of interest
        """
        indexes = self.find_ion_indexes(ion)
        atomcloud = AtomCloud()
        atomcloud.load_atomcloud(pos=self.pos[indexes,:],mass=self.mass[indexes])
        return atomcloud
    
    def filter_ions(self,ions: List[str]) -> AtomCloud:
        idx_lst = []
        for ion in ions:
            idx_lst.append(self.find_ion_indexes(ion))
        indexes = np.concatenate(idx_lst,axis=0)
        atomcloud = AtomCloud()
        atomcloud.load_atomcloud(pos=self.pos[indexes,:],mass=self.mass[indexes])
        return atomcloud

    def filter_ion_for_random_distribution(self,ion:str) -> AtomCloud:
        indexes = self.find_ion_indexes(ion)
        atomcloud = AtomCloud()
        atomcloud.load_atomcloud(pos=self.pos[indexes,:],mass=self.mass_random[indexes])
        return atomcloud


""" def _load(self,pickel):
        if pickel:
            # try load pickel file
            if os.path.isfile(self.pickelfile):
                print(f"Importing cached APT data : {self.pickelfile}")
                print(f"If you have two set of data with the same name and pickel path, please remove old pickel file to read new data, or set cache to be False.")
                self.data = pd.read_pickle(self.pickelfile)
            else:
                print(f"No cached data with name {self.name} found in directory: {self.pickel_path}. Loading data instead.")
                if os.path.isfile(self.file):
                    self.data = pd.read_csv(self.file, sep="\t",header=None)
                    self._pickel()
                else:
                    raise ValueError("Please provide correct file name!")
        else:
            if os.path.isfile(self.file):
                self.data = pd.read_csv(self.file, sep="\t",header=None)
            else:
                raise ValueError("Please provide correct file name!") """



if __name__ == '__main__':    
    pos_name = './APTdata/43334-1.6-7.POS'
    rng_name = './APTdata/43334.RRNG'
    txt_name = './APTdata/43334-1.6-5.4_Mg.txt'
    test_pos = APTData.from_pos(pos_name)
    elements, test_rng = APTData.from_rrng(rng_name)
    