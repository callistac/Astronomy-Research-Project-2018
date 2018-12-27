#following code predicts the TESS transit uncertainty in July 2019 using simulated transit techniques
#the predicted transit uncertainty for TESS from this file is used in the likelihood function in KOI_142_Analysis.py
#code makes use of Batman (Kreidberg, L. 2015)

import batman 
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import newton_krylov
g_value = 0.000295994511
M_star = 0.956

#initializing arrays 
#sampling a data point every 30 mins (but stacking transits since KOI-142 will be observed 3 times in 30 mins due to its orbital period)
thirty_min_x = np.arange(-0.2, 0.2, 0.02083/3) 
chi2_array = []
t0_array = np.linspace(-0.1, 0.1, 3000) 
new_y_array = []
t0_vs_chi2 = np.zeros(shape=(len(t0_array), 2))
thirty_min_xy_lst = []
#used for storing information from initial run of batman
flux_uncert_array = []

#creating a model lightcurve (using most probable parameters) and creating noisy "fake" data points
a_axis_au = np.cbrt((10.91657654631031)**2 *g_value*M_star/(4*np.pi**2)) 
a_axis_stellar_r = a_axis_au*214.93946938362/0.961

bat_params = batman.TransitParams()
bat_params.t0 = 0.
bat_params.per = 10.91657654631031
bat_params.rp = 0.0403
bat_params.a = a_axis_stellar_r
bat_params.inc = 88.42792792493891
bat_params.ecc = 0.05565288586927001
bat_params.w = 179.04619342089515
bat_params.u = [0.2075, 0.3785]
bat_params.limb_dark = "quadratic"

t_samples = np.linspace(-0.2, 0.2, 1300)

b_model = batman.TransitModel(bat_params, thirty_min_x, supersample_factor=7, exp_time=0.02083)
flux = b_model.light_curve(bat_params)

#finding chi squared between "observed" flux and calculated flux
np.random.seed(6)
flux_uncert = np.random.normal(0, np.sqrt(2)*1500e-6, len(thirty_min_x)) #1500 comes from Sullivan et al. paper
flux_uncert_array.append(flux_uncert)
new_y = flux + flux_uncert
new_y_array.append(new_y)
residuals = (flux - new_y)**2
chi2 = residuals / (np.sqrt(2)*1500e-6)**2 
chi2_array.append(chi2)
    
plt.plot(thirty_min_x, new_y, 'g.')
        
print(0.5*np.sum(chi2_array))
print(np.exp(-0.5*np.sum(chi2_array)))

#reinitalizing arrays to be empty
thirty_min_xy_lst = []
chi2_array = []

plt.plot(thirty_min_x, flux)
plt.xlabel("Time from central transit (days)")
plt.ylabel("Relative flux")
plt.show()

'''
altering the t0 parameter in batman to alter transit times
finds chi squared value between the previous run of batman and the new calculated flux for each specific t0
outputs a grid of t0 vs chi2 values
'''
for time in range(len(t0_array)):
    print("time", t0_array[time])
    a_axis_au = np.cbrt((10.91657654631031)**2 *g_value*M_star/(4*np.pi**2)) 
    a_axis_stellar_r = a_axis_au*214.93946938362/0.961

    bat_params = batman.TransitParams()
    bat_params.t0 = t0_array[time]
    bat_params.per = 10.91657654631031
    bat_params.rp = 0.0403
    bat_params.a = a_axis_stellar_r
    bat_params.inc = 88.42792792493891
    bat_params.ecc = 0.05565288586927001
    bat_params.w = 179.04619342089515
    bat_params.u = [0.2075, 0.3785]
    bat_params.limb_dark = "quadratic"

    t_samples = np.linspace(-0.1, 0.1, 3000) 

    b_model = batman.TransitModel(bat_params, thirty_min_x, supersample_factor=7, exp_time=0.02083) 
    flux = b_model.light_curve(bat_params)

    residuals = (flux - new_y_array)**2
    chi2 = residuals / (np.sqrt(2)*1500e-6)**2  
    chi2_array.append(chi2)
    plt.plot(thirty_min_x, new_y_array[0], 'g.')
            
    t0_vs_chi2[time, 0] = t0_array[time]
    t0_vs_chi2[time, 1] = np.sum(chi2_array)

    #reinitalizing arrays to be empty
    thirty_min_xy_lst = []
    chi2_array = []
    
 '''
finds the minimum chi2 value, subtracts that from the entire column to find chi2 values relative to the best chi2 value
converts those relative chi2 values to likelihoods (-0.5 and exponentiate)
normalizes those likelihood values
 '''
likelihood_array = np.zeros(shape=(len(t0_vs_chi2)))
transit_time = np.zeros(shape=(len(t0_vs_chi2)))
normalized_likelihood = np.zeros(shape=(len(t0_vs_chi2)))

for k in range (len(t0_vs_chi2)):
    #subtracting off the best value for chi2 so parabola has a min at zero and turning it into a likelihood
    likelihood_array[k] = (np.exp(-0.5*(t0_vs_chi2[k][1] - (np.min(t0_vs_chi2, axis=0))[1])))
    transit_time[k] = t0_vs_chi2[k][0]

#finding the norm of the likelihoods 
likelihood_array /= np.sum(likelihood_array)

'''
finding the cumulative sum of the likelihood array
subtracting off the first 16th and 84th percentile
finding the transit time values where the y value (likelihoods) crosses 0
the standard deviation is half of t_84 - t_16
'''
neg_val1 = 0.0
neg_val2 = 0.0

#finding the cumulative sum of the likelihood array
final_likelihood = np.cumsum(likelihood_array)

first_dev = final_likelihood - 0.158 
second_dev = final_likelihood - 0.842

first_func = interp1d(transit_time, first_dev)
plt.plot(first_func.x, first_func.y, '-')
plt.show()

second_func = interp1d(transit_time, second_dev)

for index1 in range(len(first_dev)):
    if first_dev[index1] < 0:
        neg_val1 = first_dev[index1]
    else:
        break
        
y1 = neg_val1
y2 = first_dev[index1]
x1 = transit_time[index1-1]
x2 = transit_time[index1]

slope1 = (y1 - y2) / (x1 - x2)
t_16 = -1*((y1 - 0) / slope1) + x1 
print("t_16", t_16)

for index2 in range(len(second_dev)):
    if second_dev[index2] < 0:
        neg_val2 = second_dev[index2]
    else:
        break
        
y84_1 = neg_val2
y84_2 = second_dev[index2]
x84_1 = transit_time[index2-1]
x84_2 = transit_time[index2]

slope2 = (y84_1 - y84_2) / (x84_1 - x84_2)
t_84 = -1*((y84_1 - 0) / slope2) + x84_1
print("t_84", t_84)

std_dev_final = 0.5*(t_84 - t_16)
print(std_dev_final)
