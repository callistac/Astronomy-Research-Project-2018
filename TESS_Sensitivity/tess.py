
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
    
'''
Function that calculates signal-to-noise ratio and inputs those SNR to our dataframe 
Input: 
    num_camp - the number of campaigns that TESS will observe the Kepler field (1=27.4 days, 2=54.8 days)
Output: 
    none, adds SNR to df
'''
def Calc_SNR(num_camp):
    for index, row in df.iterrows():
        signal = row['new_transit_depth']
        noise = row['tot_noise']
        tdur_hr = row['koi_duration']
        texp = 30 #minutes
        period = row['koi_period']
        tdur = tdur_hr * 60 
        
        #calculating the number of transits observed in a given period and number of campaigns
        num_transits = 27.4*num_camp / period
        N = int(num_transits)
        
        #the probability of detecting N+1 transits
        prob_nplus1 = num_transits%1
        
        df.at[index, 'prob_nplus1_'+ str(num_camp)] = prob_nplus1
        
        #SNR1 is the sig-noise ratio for detecting N transits
        SNR1 = (signal / noise) * np.sqrt(tdur / texp) * np.sqrt(N)
        df.at[index, 'SNR1_'+ str(num_camp)] = SNR1
        
        #SNR2 is the sig-noise tratio for detecting N+1 transits
        SNR2 = (signal / noise) * np.sqrt(tdur / texp) * np.sqrt(N+1)
        df.at[index, 'SNR2_'+ str(num_camp)] = SNR2
        
#calculates SNR for 1 or 2 campaigns (or specified value of campaigns (sectors), see Table 1 in Christ et al.
for i in range(1, 3):
    Calc_SNR(i)

#function for converting noise to a probability of detection taken from Christiansen et al. 
def Phi(x,mu=0,sigma=1):
    t = erf((x-mu)/(sigma*sqrt(2)))
    return 0.5 + 0.5*t

#taking SNR for N and N+1 for both campaigns and converting them to probabilities of detection
#tess has a 7.1 signal detection threshold, but best case scenario would be a 3 sigma signal detection threshold
for index, row in df.iterrows():
    snratio1_camp1 = row['SNR1_1']
    snratio2_camp1 = row['SNR2_1']
    snratio1_camp2 = row['SNR1_2']
    snratio2_camp2 = row['SNR2_2']
    prob_num_transits1 = row['prob_nplus1_1']
    prob_num_transits2 = row['prob_nplus1_2']
    
    #this line was used to find total SNR for harvard people
    df.at[index, 'total_SNR1'] = (snratio1_camp1*(1 - prob_num_transits1)) + (snratio2_camp1*(prob_num_transits1))
    df.at[index, 'total_SNR2'] = (snratio1_camp2*(1 - prob_num_transits2)) + (snratio2_camp2*(prob_num_transits2))
    
    period = row['koi_period']
    prob1_camp1 = Phi(snratio1_camp1, mu=3.0, sigma=1.0) #mu = signal detection threshold, using either 3 or 7.1
    prob2_camp1 = Phi(snratio2_camp1, mu=3.0, sigma=1.0)
    
    prob1_camp2 = Phi(snratio1_camp2, mu=3.0, sigma=1.0)
    prob2_camp2 = Phi(snratio2_camp2, mu=3.0, sigma=1.0)
    
    #tot_prob = probability of N transit* probability of detecting one transit + probability of N+1 transits* prob of detecting N+1 transits
    tot_prob_camp1 = prob1_camp1*(1 - prob_num_transits1) + prob2_camp1*(prob_num_transits1)
    tot_prob_camp2 = prob1_camp2*(1 - prob_num_transits2) + prob2_camp2*(prob_num_transits2)
    
    df.at[index, 'tot_prob_camp1'] = tot_prob_camp1
    df.at[index, 'tot_prob_camp2'] = tot_prob_camp2
 
#finding out how many of these prob vals are nan
nan_df = df['tot_prob_camp1'][~np.isfinite(df['tot_prob_camp1'])]
print(len(nan_df))

#finding any planets that have either a missing tess magnitude or depth
missing_tess_mag = []
missing_depth = []
for index, row in df.iterrows():
    if index in list(nan_df.index.values):
        if (np.isnan(df.loc[index]['koi_depth'])):
            missing_depth.append(df.loc[index]['kepid'])
        if (np.isnan(df.loc[index]['tess_magnitude'])):
            missing_tess_mag.append(df.loc[index]['kepid'])

print("missing tess mag", missing_tess_mag)
print("missing depth", missing_depth)
for kk in range(len(missing_tess_mag)):
    if missing_tess_mag[kk] in missing_depth:
        print(missing_tess_mag[kk])
        
#if a planet is missing a depth, then get rid of it
for index, row in df.iterrows():
    if row['kepid'] in missing_depth:
        df = df.drop(labels=index)
        
