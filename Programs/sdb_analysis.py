#! /usr/bin/env python
### Script containing varous plotting functions for splitting Measurements
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import cartopy.crs as cart
#import cartopy
import matplotlib.gridspec as gridspec
import obspy
from obspy import taup



class SDB:
    """
    Class to hold a Splitting Database (and Pierce Points) and generate a suite of useful plots based off the data.
    """
    def __init__(self,sdb,Q_threshold=None,gcarc_threshold=None):

        self._raw = pd.read_csv('{}.sdb'.format(sdb),delim_whitespace=True)
        #Load raw data from provided sdb file. This is going to be a hidden file as I will parse out useful columns to new attributes depending of provided kwargs

         #if kwargs are none:
        self.data = self._raw
