import numpy as np
import pandas as pd

from aptca.core.APTData import APTData
from hdbscan import HDBSCAN


class Clusters:
    """
    Container for cluster information: atom position, persistence, probability, etc.
    """
    def __init__(self, aptData: APTData, pos:np.ndarray, clusterer:HDBSCAN):
        raise NotImplementedError()





    

    

    
    
