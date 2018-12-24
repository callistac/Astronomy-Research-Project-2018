# If you make use of this code, please cite Christ, C. N., Montet, B. T., & Fabrycky, D. C. 2018, arXiv:1810.02826

# The following code is used to determine how much TESS will improve our measurements of KOI-142 (Kepler-88)

import numpy as np
import matplotlib.pyplot as pyplot
import scipy.optimize as op
import os
import emcee
import corner

ndim, nwalkers = 13, 150 
g_value = 0.000295994511
M_star = 0.956

G_guess = 0.000295994511
Mstar_guess = 0.956
M1_guess = 4.34568221e-05
P1_guess = 10.915996767 
E1_guess = 0.056137334 
i1_guess = 89.116000000
LongNode1_guess = 0.000000000 
W1_guess = -179.794778327+360
mean1_guess = 261.169991641
M2_guess = 6.113274916e-04
P2_guess = 22.267768543 
E2_guess = 0.056964888 
i2_guess = 85.666179336
LongNode2_guess = -0.583597905+360 
W2_guess = 0.081068179 
mean2_guess = 335.896661714

#initializing arrays
tot_transit = np.zeros((355, 10000))
mean_transit = np.zeros(355)
stdeviation_transit = np.zeros(355)
chi_squared = np.zeros((1, 120))
miss_epoch = np.array([4, 10, 21, 24, 55, 56, 86, 91, 95, 107, 108, 118, 120, 124]) - 1

#loading KOI-142 Kepler data
txt_file = np.loadtxt("/Users/Callista/Documents/Github/TTVFast/c_version/142_TTVs.txt")
TOT = txt_file[:, 1]
obs_times = txt_file[:, 2] + txt_file[:, 3] / 1440.0
error =  txt_file[:,4] / 1440.0

#initializing variables to solve for chi squared in transit duration
R_star_AU = 0.961*(1/214.9394693836)
delta = 0.039
obs_times_dur = 3.28 
error_dur = 0.185

'''
Goal is to find the posterior distribution (and best fit values) of the parameters for KOI-142
*Note* Same process was done in Nesvorny et al. (2013), but they did not post their posteriors and so we redid their analysis with all Kepler data

Creating likelihood function that will be used for MCMC analysis (chi squared between observed and predicted values by TTVFast)
  input: theta (parameter values we are looking to alter) additional unchanging parameters (G_VAL and Mstar)
  output: function will return a chi squared value for a given set of initial parameters (measures how well the observed distribution of data fits with the distribution that is expected if the variables are independent)-goodness of fit
'''
def likelihood(theta, G_VAL, Mstar):
    Mpb, Periodb, Eb, Ib, Wb, Meanb, Mpc, Periodc, Ec, Ic, Wc, Longnodec, Meanc = theta
    Longnodeb = 0.0
    #creating/writing .in file
    filename = ("/Users/Callista/Documents/GitHub/infiles/TTVs0" + ".in")
    infile = open(filename, 'w')
    infile.write("%.11f\n%.11f\n%.11f\n" % (G_VAL, Mstar, Mpb))
    infile.write("%.11f %.11f %.11f %.11f %.11f %.11f\n" % (Periodb, Eb, Ib, Longnodeb, Wb, Meanb))
    infile.write("%.11f\n" % Mpc)
    infile.write("%.11f %.11f %.11f %.11f %.11f %.11f\n" % (Periodc, Ec, Ic, Longnodec, Wc, Meanc))
    infile.close()

    #creating/writing setup file
    setupfilename = ("/Users/Callista/Documents/Github/setupfiles/new_setupfile0")
    new_setupfile = open(setupfilename, 'w')
    new_setupfile.write("%s\n %.8f\n %.3f\n %d\n %d\n %d\n" % (filename, 54.675215, 0.54, 3950, 2, 0))
    new_setupfile.close()
    os.system(("./run_TTVFast" + " " + setupfilename + " /Users/Callista/Documents/Github/KOI142_files/final_files0" + " RV_file RV_out"))

    tmp_array = np.loadtxt("/Users/Callista/Documents/Github/KOI142_files/final_files0")
    planet = tmp_array[:,0]
    epoch = tmp_array[:,1]
    time = tmp_array[:,2]
    Vsky = tmp_array[:, 4]

    planet_1 = planet[planet == 0]
    epoch_1 = epoch[planet == 0]
    time_1 = time[planet == 0]
    Vsky_1 = Vsky[planet == 0]
    time_1 += 0./1440

    #finding chi-squared values
    if time_1[0] < 60:
        slice_time = time_1[1:]
        epoch_1 = np.delete(epoch_1, miss_epoch)
        time_1 = np.delete(slice_time, miss_epoch)
        chi_squared[i, :] = ((time_1[:120] - obs_times[1:])/ (error[1:]))**2
        chi2_tess = ((time_1[-12] - 3823.0106513) / 0.01174660432206573)**2 
    else:
        epoch_1 = np.delete(epoch_1, miss_epoch)
        time_1 = np.delete(time_1, miss_epoch)
        chi_squared[i, :] = ((time_1[:120] - obs_times[1:]) / (error[1:]))**2
        chi2_tess = ((time_1[-12] - 3823.0106513) / 0.01174660432206573)**2 

    #conversions
    Ib_rad = Ib * 0.0174533 
    Wb_rad = Wb * 0.0174533
    
    #finding transit duration vals 
    semi_axis = np.cbrt(Periodb**2 * g_value*M_star/(4*np.pi**2))
    b_val = ((semi_axis**2) * (np.cos(Ib_rad))**2 * (1 / R_star_AU**2) * ((1 - Eb**2)/(1 + Eb * np.sin(Wb_rad)))**2)
    
    if b_val > 1:
        print(np.inf)
        return np.inf
    
    Tdur = ((2 * (1 + delta) * R_star_AU / Vsky_1) * np.sqrt(1 - b_val))

    #converting transit duration values to hours and finding chi squared val
    Tdur = Tdur[:1]
    Tdur_hours = Tdur * (24.0)
    chi2_val = ((Tdur_hours - obs_times_dur) / (error_dur))**2
      
    #summing chi_vals
    sum_matrix = np.sum(chi_squared, axis=1)
    ind = np.unravel_index(np.argmin(sum_matrix, axis=None), sum_matrix.shape)
    print(np.sum(0.5*(sum_matrix[ind] + chi2_val + chi2_tess))) 
    return 0.5*(sum_matrix[ind] + chi2_val + chi2_tess) 

#make sure chi squared values significantly are reduce by end of optimization
bnds = ((0, 0.006), (10.90, 10.93), (0, 1), (60, 90), (0, 360), (0, 360), (0, 0.006), (22.20, 22.35), (0, 1), (60, 120), (180, 540), (-180, 180), (0, 360))
result = op.minimize(likelihood, [M1_guess, P1_guess, E1_guess, i1_guess,
                          W1_guess, mean1_guess, M2_guess, P2_guess, E2_guess, i2_guess, LongNode2_guess,
                          W2_guess, mean2_guess], args=(g_value, M_star), method="L-BFGS-B", bounds= bnds)

M1, P1, E1, I1, W1, Mean1, M2, P2, E2, I2, Longnode2, W2, Mean2 = result["x"] 

'''
defining a prior for the total probability function
'''
def prior(theta):
    Mpb, Periodb, Eb, Ib, Wb, Meanb, Mpc, Periodc, Ec, Ic, Wc, Longnodec, Meanc = theta
    #force range to be above or below 90 for inclination
    if ((0 < Mpb < 0.006) and (10.90 < Periodb < 10.93) and (0 < Eb < 1) and 
        (60 < Ib < 90) and (0 < Wb < 360) and (0 < Meanb < 360) and (0 < Mpc < 0.006) and 
        (22.20 < Periodc < 22.35) and (0 < Ec < 1) and (60 < Ic < 120) and (180 < Wc < 540) and (-180 < Longnodec < 180)
        and (0 < Meanc < 360)): 
            return 0.0
    else:
        return -np.inf

'''
Defining a total probability function 
'''
def tot_prob(theta, G_VAL, Mstar):
    lp = prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    else:
        return lp + -1*likelihood(theta, G_VAL, Mstar)
 
