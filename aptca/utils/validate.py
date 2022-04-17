from __future__ import annotations

import os
from os.path import abspath, dirname, join
import numpy as n

def file_exists(fpath: str) -> str:
    """
    Validate that a file exists

    :param fpath: the file path to validate existence
    """
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"The path {fpath} does not exist")
    elif not os.path.isfile(fpath):
        raise IOError(f"The path {fpath} is not a file")

    return fpath