#the following code is used to create figures to visualize how TESS in conjunction with Kepler improves our knowledge of KOI-142
#see KOI_142_Analysis.py file for important chains used in the following code

import matplotlib.pyplot as plt
import numpy as np

#initializing plot layout
plt.rcParams['axes.linewidth']=3
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.minor.width'] = 2
plt.rcParams['ytick.minor.width'] = 2
plt.rc('xtick.major', size=8, pad=8)
plt.rc('xtick.minor', size=6, pad=5)
plt.rc('ytick.major', size=8, pad=8)
plt.rc('ytick.minor', size=6, pad=5)

#chains from KOI_142_analysis.py
c = np.load('sampler_chains2.npy')
tess_chains = np.load('sampler_chains2_tess.npy')

#creating an array of dates starting from 2009 until 2028
time_array28 = np.zeros(shape=(641,))
for time28 in range(len(time_array28)):
    if (time28==0):
        time_array28[time28] = 2009.25
    else:
        time_array28[time28] = time_array28[time28-1] + 10.915996767 / 365

#initializing arrays
chain_transits28 = np.zeros(shape=(1500000, 641)) 
tess_chain_transits28 = np.zeros(shape=(1500000, 641))

'''
Calculates transit times 
Input:
    full_chain (array of chains from MCMC analysis)
    specific_chains (tuple, which chains to include in analysis out of full_chain)
    ending_date (int, which date to end analysis, in this case our ending date is 2028)
Output:
    transit times
'''
def calc_transit_times(full_chain, specific_chains, ending_date):
    params = full_chain[specific_chains[0], specific_chains[1], :]
    filename = ("/Users/Callista/Documents/GitHub/infiles2/TTVs0" + ".in")
    infile = open(filename, 'w')
    infile.write("%.11f\n%.11f\n%.11f\n" % (g_value, M_star, params[0]))
    infile.write("%.11f %.11f %.11f %.11f %.11f %.11f\n" % (params[1], params[2], params[3], 0.0, params[4], params[5]))
    infile.write("%.11f\n" % params[6])
    infile.write("%.11f %.11f %.11f %.11f %.11f %.11f\n" % (params[7], params[8], params[9], params[11], params[10], params[12]))
    infile.close()
    
    #creating/writing setup file
    setupfilename = ("/Users/Callista/Documents/Github/setupfiles2/new_setupfile0")
    new_setupfile = open(setupfilename, 'w')
    new_setupfile.write("%s\n %.8f\n %.3f\n %d\n %d\n %d\n" % (filename, 54.675215, 0.54, ending_date, 2, 0))
    new_setupfile.close()
    os.system(("./run_TTVFast" + " " + setupfilename + " /Users/Callista/Documents/Github/KOI142_files2/final_files0" + " RV_file RV_out"))
    
    transit_times_file = np.loadtxt("/Users/Callista/Documents/Github/KOI142_files2/final_files0")
    
    planet = transit_times_file[:,0]
    epoch = transit_times_file[:,1]
    time = transit_times_file[:,2]
        
    planet_1 = planet[planet == 0]
    epoch_1 = epoch[planet == 0]
    time_1 = time[planet == 0]
    time_1 += 0./1440
    return time_1

#stepping through each set of parameters in the array of chains and calculating transit times until 2028 for each set of parameters
#takes all of the possible sets of transits based on the given sets of parameters and inputs these times into a new array
ind_time = 0 
for i in range(len(c)):
    for j in range(len(c[0])):
        ttr = calc_transit_times(c, (i, j), 7067)
        tess_transits = calc_transit_times(tess_chains, (i, j), 7067)
        chain_transits28[ind_time] = ttr
        tess_chain_transits28[ind_time] = tess_transits   
        ind_time+=1

#finds the standard deviation of the transit times (column of array)
#given a specific time, std_array2 and std_array will be our uncertainity in the transit times with and without TESS, respectively
std_array = np.std(chain_transits28, axis=0)*1440
std_array2 = np.std(tess_chain_transits28, axis=0)*1440

#plot figure
plt.figure(figsize=(18,7))
plt.gcf().subplots_adjust(left=0.17, bottom=0.17, right=0.94, top=0.94, wspace=0.0, hspace=0.0)

plt.axvline(x=2019.54166667, linewidth=2.5, color='r', alpha=0.8)
plt.axvline(x=2019.5, linewidth=0.5, color='k', alpha=0.8)
plt.axvline(x=2019.583333333, linewidth=0.5, color='k', alpha=0.8)
plt.plot(time_array28, std_array, 'g.', time_array28, std_array2, 'b.', markersize=9)

plt.gcf().subplots_adjust(left=0.17, bottom=0.17, right=0.94, top=0.94, wspace=0.0, hspace=0.0)
plt.gca().set_ylim(bottom=0)
plt.gca().set_xlim(right=np.max(time_array28))

plt.xlabel('Time (years)', fontsize=24)
plt.ylabel('Uncertainty in transit time (minutes)', fontsize=20)
plt.xticks(fontsize=19)
plt.yticks(fontsize=19)
plt.savefig('final_transit_uncert_new.pdf', bbox_inches='tight')
