
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

#falsely labeled planet
pos1 = np.where(df['kepid'] == 8110757)
df.at[6993, 'koi_disposition'] = 'FALSE POSITIVE'

#only using confirmed or candidate systems in the following analysis
df = df.loc[df['koi_disposition'].isin(['CONFIRMED', 'CANDIDATE'])]

#file of kepler contam ratios
file1 = open('kepler_fov_search.txt')
df2 = pd.read_csv(file1, names=[ 'Kepler_ID', 'Crowding_season_0', 'Crowding_season_1', 'Crowding_season_2', 'Crowding_season_3'], skiprows=2)
file1.close()

'''
the following code is all pre-processing of the TESS/Kepler data to ensure that we are not missing any information for each system
'''
#recording systems missing from the kepler contam file (should be none)
#taking median value of contam ratio and inputting in dataframe for kepler
missing_systems = []
for index, row in df.iterrows():
    idnum = row['kepid']
    loc = np.where(df2['Kepler_ID'] == idnum)
    if (len(loc[0]) == 0):
        missing_systems.append(idnum)
    else:
        rowdf2 = df2.iloc[loc[0]]
        crowding_array = np.array([rowdf2.iloc[0]['Crowding_season_0'], rowdf2.iloc[0]['Crowding_season_1'], rowdf2.iloc[0]['Crowding_season_2'], rowdf2.iloc[0]['Crowding_season_3']])
        contam_val = np.nanmedian(crowding_array)
        df.at[index, 'kep_contam_ratio'] = contam_val 
        
#tess contam ratio files 
df3 = pd.read_csv('MAST_Crossmatch_TIC.csv', skiprows=4)
df4 = pd.read_csv('MAST_Crossmatch_TIC_missing_Tmag.csv', skiprows=4)

#recording missing systems from tess contam file (should be a fair amount) 
#search for those missing systems in another dataset and input the tess magnitude/contamratio into df
#if not missing, input tess mag/contamratio into df
#should have no missing_systems2 after running this code
missing_systems2 = []
nan_systems = []
for index, row in df.iterrows():
    idnum = row['kepid']
    loc = np.where(df3['Kepler_ID'] == idnum)
    if (len(loc[0]) == 0):
        loc1 = np.where(df4['Kepler_ID'] == idnum)
        if (len(loc1) == 0):
            missing_systems2.append(idnum)
        else: 
            row2df = df4.iloc[(loc1[0][0])]
            missing_tmag = row2df['Tmag']
            missing_contamratio = row2df['contratio']
            if np.isnan(missing_contamratio):
                nan_systems.append(idnum)
                df.at[index, 'tess_magnitude'] = missing_tmag
                df.at[index, 'tess_contam_ratio'] = 0.0
            else:
                df.at[index, 'tess_magnitude'] = missing_tmag
                df.at[index, 'tess_contam_ratio'] = missing_contamratio
    else:
        rowdf2 = df3.iloc[loc[0]]
    
        contam_val = rowdf2.iloc[0]['contratio']
        tess_mag = rowdf2.iloc[0]['Tmag']
        
        if (np.isnan(contam_val) == True):
            nan_systems.append(idnum)
            contam_val = 0.0
            
        df.at[index, 'tess_contam_ratio'] = contam_val
        df.at[index, 'tess_magnitude'] = tess_mag
        

#contaminating transit depts with the tess contam ratio (first have to convert t_contam_ratio to the same form as k_contam_ratio)
#input new transit depth into df
for index, row in df.iterrows():
    koidepth = row['koi_depth'] 
    k_contam_ratio = row['kep_contam_ratio']
    t_contam_ratio = row['tess_contam_ratio']
    t_contam_ratio = 1 / (1 + t_contam_ratio)
    koidepth = koidepth * t_contam_ratio  
    df.at[index, 'new_transit_depth'] = koidepth
    
 #coming up with the equation that will convert a tess magnitude to total noise that tess will observe
#used sullivan et al's plot to come up with the equation for total noise

for index, row in df.iterrows():
    mag = row['tess_magnitude']
    
    star_noise_log = 0.2235*mag + 0.0565
    sky_noise_log = 0.3347*mag - 1.6776
    read_noise_log = 0.3486*mag - 1.940
    syst_noise_log = np.log10(60)
    
    star_noise_r = 10**star_noise_log
    sky_noise_gn = 10**sky_noise_log
    read_noise_gy = 10**read_noise_log
    syst_noise_b = 10**syst_noise_log
    
    tot_noise = np.sqrt(star_noise_r**2 + sky_noise_gn**2 + read_noise_gy**2 + syst_noise_b**2)
    tot_noise *= np.sqrt(2)
    df.at[index, 'tot_noise'] = tot_noise 
