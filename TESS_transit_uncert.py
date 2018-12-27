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
