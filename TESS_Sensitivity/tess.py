
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from scipy.special import erf

#reading in data from the NASA Exoplanet Archive
df = pd.read_csv('exoplanet_archive_complete_version4.csv', skiprows=53)

#including only the following columns/data:
df = df.filter(items=['rowid', 'kepid', 'kepoi_name', 'kepler_name', 'koi_disposition', 'koi_depth', 'koi_duration', 'koi_prad', 'koi_period', 'koi_period_err1', 'koi_time0bk', 'koi_time0bk_err1', 'koi_smass', 'koi_srad', 'koi_sma', 'koi_eccen', 'koi_incl', 'koi_longp'])

#adding in columns that we will later have data for
df['kep_contam_ratio'] = pd.Series(np.nan, index=df.index)
df['tess_contam_ratio'] = pd.Series(np.nan, index=df.index)
df['new_transit_depth'] = pd.Series(np.nan, index=df.index)
df['tess_magnitude'] = pd.Series(np.nan, index=df.index)
df['tot_noise'] = pd.Series(np.nan, index=df.index)
df['SNR1_1'] = pd.Series(np.nan, index=df.index)
df['SNR2_1'] = pd.Series(np.nan, index=df.index)
df['SNR1_2'] = pd.Series(np.nan, index=df.index)
df['SNR2_2'] = pd.Series(np.nan, index=df.index)

