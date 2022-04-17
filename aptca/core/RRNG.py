"""
Map elements (ions) to mass-to-charge ratio ranges
"""

from distutils.log import warn
import os
import re
from typing import Set
import numpy as np
import pandas as pd

from aptca.utils import validate
from aptca.utils.logging import log


class RRNG:

    def __init__(self,
                elements : np.ndarray = None,
                rrng     : np.ndarray = None):
        self._elements = elements
        self._rrng = rrng
        self._loaded = False
        if self._rrng != None:
            self._loaded = True
            self._ions = self.get_unique_ions(self.rrng)
    
    def load_rrng_fromfile(self,rrng_file:str):
        validate.file_exists(rrng_file)
        self._elements, self._rrng = RRNG.from_rrng(rrng_file)
        self._ions = self.get_unique_ions(self._rrng)
        self._loaded = True

    def get_unique_ions(self,rrng):
        ions = set()
        for rng in rrng:
            ions.add(rng[3])
        return np.sort(np.array(list(ions),dtype=str))    

    def get_range_for_ion(self,ion:str) -> np.ndarray:
        ranges = []
        if self._loaded:
            for rng in self._rrng:
                if rng[3] == ion:
                    ranges.append([rng[0],rng[1]])
            return np.array(ranges)
        else:
            log.error("Please load RRNG file first!")
            raise ValueError("RRNG object empty!")
    
    def get_ion_from_mass(self,mass: float) -> str:
        for rng in self._rrng:
            if mass >= rng[0] and mass <= rng[1]:
                return rng[3]
        else:
            return '-1'

    @property
    def elements(self) -> np.ndarray:
        return self._elements

    @property
    def rrng(self) -> np.ndarray:
        return self._rrng
    
    @property
    def ions(self) -> np.ndarray:
        return self._ions

    @classmethod
    def from_rrng(cls,rrng_fname):
        """
        Loads a .rrng file (IVAS format). Returns an array contains unique 'ions', and an array for 'ranges'.
        Range file format for IVAS is inconsistent among known and unknown range names.

        For known range name (contains only ion names on element table), format is:

            Range1=low_m high_m Vol:some_float ion_species:number (ion_species_2:number ...) Color:hexdecimal_num

        For unknown range names (any name given that not on element table):
            Range1=low_m high_m Vol:some_float Name:some_name Color:hexdecimal_num
        
        rrng_fname:
            filename for the range file.
        
        return:
            (ions, rrng): ions is a numpy array for all ion species in the range file; 
            rrng is a structured numpy array [(low_m, high_m, vol, ion_type,.., color), ...]
            dt = np.dtype([('range_low', '>f4'), ('range_high', '>f4'), ('vol', '>f4'), ('ion_type', 'U16'), ('color', 'U16')])
        """
        log.info("Reading RRNG file: {}".format(os.path.basename(rrng_fname)))
        with open(rrng_fname, 'r') as f:
            rf = f.readlines()
        
        # pattern is a compiled regular experssion pattern to match range file format. works for both known unknown names.
        # pattern maching for strings: 
        #       r: raw string literal flag, means '\' won't be treated as escape character.
        #       'Ion([0-9]+)=([A-Za-z0-9]+)': it matches the ion types in first section such as ion1=C, ion2=Ti1O1,...
        #                                     'Ion' will match exactly character; (...) means a group, which will be retrived; [0-9] means a set of characters, in
        #                                     this case, any numeric character between 0 to 9; + means one or more repetition of the preceding regular experssion.
        #       '|' means 'or', that is either the precede one or the trailing one is matched.
        #       '-?': to match 0 or more minus sign, just put here for test purpose (some times I use -1 as M/Z ratio for noise points in test data)
        #       '\d+.\d': to match a float number, such as 123.456
        #       '([a-zA-Z0-9:a-zA-Z0-9 ]+)' : to match ion types, such as 'Ti:1 O:1' or 'Name:Cluster1'. Note the last space within square paranthesis is important.
        #       'Color:([A-Za-z0-9]{6})': to match color hexdecimal number. It matches exactly 6 characters.
        pattern = re.compile(r'Ion([0-9]+)=([A-Za-z0-9]+)|Range([0-9]+)=(-?\d+.\d+) (-?\d+.\d+) +Vol:(\d+.\d+) +([a-zA-Z0-9_+:a-zA-Z0-9 ]+) +Color:([A-Za-z0-9]{6})')
            
        elements = []
        rrngs = []
        for line in rf:
            match_re = pattern.search(line)
            if match_re:
                if match_re.groups()[0] is not None:
                    elements.append(list(match_re.groups())[1])
                else:
                    rrngs.append(match_re.groups()[3:])
                    
        dt = np.dtype([('range_low', '>f4'), ('range_high', '>f4'), ('vol', '>f4'), ('ion_type', 'U16'), ('color', 'U16')]) # Note that in python 3 all strings are unicode now.
        rrngs = np.array(rrngs, dtype=dt)
        elements = np.array(elements)
        
        pattern_unknown_ion = re.compile(r'(?<=Name)([A-Za-z0-9]+)') # To obtain correct ion name for unknown ions in range file. Separate the 'Name' from e.g. 'Name:Cluster1'.

        # To further process ion types, reorganize like 'Ti:1 O:1' to 'Ti1O1'.
        for idx in range(len(rrngs)):
            rrngs[idx][3] = rrngs[idx][3].replace(':', '')
            rrngs[idx][3] = rrngs[idx][3].replace(' ', '')

            n = pattern_unknown_ion.search(rrngs[idx][3])
            if n:
                rrngs[idx][3] = n.groups()[0]

        # check the validity of range, there should be no interval overlap.
        sorted_rng_idx = np.argsort(rrngs['range_low'], kind='mergesort')
        range_low = rrngs['range_low']
        range_low = range_low[sorted_rng_idx]
        range_high = rrngs['range_high']
        range_high = range_high[sorted_rng_idx]
        assert np.all(range_low < range_high), 'Invalid range file: range overlap detected.'

        return elements, rrngs[sorted_rng_idx] 