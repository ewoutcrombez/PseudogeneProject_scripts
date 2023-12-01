import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import lmfit

y_filt_psg = np.array([1195, 811, 409, 883, 914, 911, 28, 37, 82, 13, 18]) # Number of pseudogenes (filtered Ks)
y_filt_ap = np.array([17865,8783,2806,15756,13620,9269,2380,1606,843,800,811])/20 # Number of anchorpairs (filtered Ks)
y_psg = np.array([2248, 1339, 750, 1551, 1640, 1591, 109, 241, 263, 86, 179]) # Number of pseudogenes (unfiltered Ks)
y_ap = np.array([56194, 13604, 4786, 23388, 31785, 19942, 7512, 7643, 3475, 3641, 2061])/30 # Number of anchorpairs (unfiltered Ks)
x_timing = np.array([13.59, 18.32, 20.40, 26.78, 26.78, 34.73, 50.07, 59.56, 66.23, 69.05, 120]) # Timing of WGM events
x_Ks = np.array([0.1, 0.2, 0.2, 0.4, 0.4, 0.275, 0.8, 0.7, 0.8, 0.95, 1.25]) # Ks of WGM events (Vanneste)
y_filt_ap_without_bra_and_sbi = np.array([17865,8783,2806,13620,9269,2380,1606,843,811])/20 # Number of anchorpairs (filtered Ks)
y_filt_psg_without_bra_and_sbi = np.array([1195, 811, 409, 914, 911, 28, 37, 82, 18]) # Number of pseudogenes (filtered Ks)
x_timing_without_bra_and_sbi = np.array([13.59, 18.32, 20.40, 26.78, 34.73, 50.07, 59.56, 66.23, 120]) # Timing of WGM events
x_Ks_without_bra_and_sbi = np.array([0.1, 0.2, 0.2, 0.4, 0.275, 0.8, 0.7, 0.8, 1.25]) # Ks of WGM events (Vanneste)

# The double exponential function that we want to fit our data to
# This corresponds to the two-phase model (see https://doi.org/10.1073/pnas.1507669112)
def func_double_exp(x, A, q, a, c):
    return A * ((1 - q) * np.exp(-a * x) + q * np.exp(-c * x))

# The exponential function that we want to fit our data to
# This corresponds to the passive or random loss model
def func_exp(x, a, b):
    return a * np.exp(b*x)

def fit_curves(x, y, name, p0_exp = [2000, 0.02], p0_double_exp = [2000, 0.5, 0.01, 0.0001]):
    # Curve fit the test data with func_exp
    fittedParameters, pcov = curve_fit(func_exp, x, y, p0 = p0_exp, maxfev = 1000000)
    modelPredictions = func_exp(x, *fittedParameters)
    ## Calculate R-squared
    rss = sum((y - modelPredictions)**2)
    ss_tot = sum((y - np.mean(y))**2)
    r_squared_exp = 1 - (rss / ss_tot)

    ## Obtain parameters from the fitted curve
    a_exp = fittedParameters[0]
    b_exp = fittedParameters[1]

    # Curve fit the test data with func_double_exp
    fittedParameters, pcov = curve_fit(func_double_exp, x, y, p0 = p0_double_exp, maxfev = 1000000)
    modelPredictions = func_double_exp(x, *fittedParameters)
    ## Calculate R-squared
    rss = sum((y - modelPredictions)**2)
    ss_tot = sum((y - np.mean(y))**2)
    r_squared_double_exp = 1 - (rss / ss_tot)
    
    ## Obtain parameters from the fitted curve
    A_double_exp = fittedParameters[0]
    q_double_exp = fittedParameters[1]
    a_double_exp = fittedParameters[2]
    c_double_exp = fittedParameters[3]

    # Also fit the data using lmfit: exponential
    model = lmfit.models.ExpressionModel("a * exp(b*x)")
    ## Set initial parameter values
    params = model.make_params(a = p0_exp[0], b = p0_exp[1])
    ## Fit the model to the data
    fit = model.fit(y, params, x=x)
    print(fit.fit_report())
    ## calculate R-squared
    rss = sum((y - fit.best_fit)**2)
    ss_tot = sum((y - np.mean(y))**2)
    r_squared_exp_lmfit = 1 - (rss / ss_tot)

    ## Obtain parameters from the fitted curve
    a_exp_lmfit = fit.best_values['a']
    b_exp_lmfit = fit.best_values['b']

    # Also fit the data using lmfit: double exponential
    model = lmfit.models.ExpressionModel("A * ((1-q)*exp(-a*x) + q*exp(-c*x))")
    ## Set initial parameter values
    params = model.make_params(A = p0_double_exp[0], q = p0_double_exp[1], a = p0_double_exp[2], c = p0_double_exp[3])
    ## constrain q to be between 0 and 1 (0 <= q <= 1)
    params['q'].min = 0
    params['q'].max = 1
    ## Fit the model to the data
    fit = model.fit(y, params, x=x)
    print(fit.fit_report())
    ## calculate R-squared
    rss = sum((y - fit.best_fit)**2)
    ss_tot = sum((y - np.mean(y))**2)
    r_squared_double_exp_lmfit = 1 - (rss / ss_tot)

    ## Obtain parameters from the fitted curve
    A_double_exp_lmfit = fit.best_values['A']
    q_double_exp_lmfit = fit.best_values['q']
    a_double_exp_lmfit = fit.best_values['a']
    c_double_exp_lmfit = fit.best_values['c']

    result_string = name + "\t" + str(r_squared_exp) + "\t" + str(a_exp) + "\t" + str(b_exp) + "\t" + str(r_squared_double_exp) + "\t" +\
        str(A_double_exp) + "\t" + str(q_double_exp) + "\t" + str(a_double_exp) + "\t" + str(c_double_exp) + "\t" +\
        str(r_squared_exp_lmfit) + "\t" + str(a_exp_lmfit) + "\t" + str(b_exp_lmfit) + "\t" +\
        str(r_squared_double_exp_lmfit) + "\t" + str(A_double_exp_lmfit) + "\t" + str(q_double_exp_lmfit) + "\t" +\
        str(a_double_exp_lmfit) + "\t" + str(c_double_exp_lmfit) + "\n"
    return result_string

string = "comparison\tr_squared_exp\ta_exp\tb_exp\tr_squared_double_exp\tA_double_exp\tq_double_exp\ta_double_exp\tc_double_exp\t" +\
    "r_squared_exp_lmfit\ta_exp_lmfit\tb_exp_lmfit\tr_squared_double_exp_lmfit\tA_double_exp_lmfit\tq_double_exp_lmfit\t" +\
    "a_double_exp_lmfit\tc_double_exp_lmfit\n"
psg_vs_Ks = fit_curves(x_Ks, y_psg, "psg_vs_Ks", p0_double_exp = [2000, 0.8, 0.1, 0.00001])
string += psg_vs_Ks
psg_vs_Ks_filtered = fit_curves(x_Ks, y_filt_psg, "psg_vs_Ks_filtered", p0_double_exp = [2000, 0.8, 0.1, 0.000001])
string += psg_vs_Ks_filtered
psg_vs_timing = fit_curves(x_timing, y_psg, "psg_vs_timing", p0_double_exp = [2000, 0.8, 0.1, 0.000001])
string += psg_vs_timing
psg_vs_timing_filtered = fit_curves(x_timing, y_filt_psg, "psg_vs_timing_filtered", p0_double_exp = [2000, 0.8, 0.1, 0.000001])
string += psg_vs_timing_filtered
APs_vs_Ks = fit_curves(x_Ks, y_ap, "APs_vs_Ks", p0_double_exp = [2000, 0.5, 0.01, 0.01])
string += APs_vs_Ks
APs_vs_Ks_filtered = fit_curves(x_Ks, y_filt_ap, "APs_vs_Ks_filtered", p0_double_exp = [2000, 0.01, 1, 0.001])
string += APs_vs_Ks_filtered
APs_vs_timing = fit_curves(x_timing, y_ap, "APs_vs_timing", p0_double_exp = [2000, 0.5, 0.01, 0.0001])
string += APs_vs_timing
APs_vs_timing_filtered = fit_curves(x_timing, y_filt_ap, "APs_vs_timing_filtered", p0_double_exp = [2000, 0.5, 0.1, 0.0001])
string += APs_vs_timing_filtered
# with Brassica rapa and Sorghum bicolor removed
APs_vs_Ks_without_bra_and_sbi = fit_curves(x_Ks_without_bra_and_sbi, y_filt_ap_without_bra_and_sbi, "APs_vs_Ks_without_bra_and_sbi", p0_double_exp = [2000, 0.01, 1, 0.00001])
string += APs_vs_Ks_without_bra_and_sbi
psg_vs_Ks_without_bra_and_sbi = fit_curves(x_Ks_without_bra_and_sbi, y_filt_psg_without_bra_and_sbi, "psg_vs_Ks_without_bra_and_sbi", p0_double_exp = [2000, 0.8, 0.1, 0.000001])
string += psg_vs_Ks_without_bra_and_sbi
APs_vs_timing_without_bra_and_sbi = fit_curves(x_timing_without_bra_and_sbi, y_filt_ap_without_bra_and_sbi, "APs_vs_timing_without_bra_and_sbi", p0_double_exp = [2000, 0.5, 0.1, 0.0001])
string += APs_vs_timing_without_bra_and_sbi
psg_vs_timing_without_bra_and_sbi = fit_curves(x_timing_without_bra_and_sbi, y_filt_psg_without_bra_and_sbi, "psg_vs_timing_without_bra_and_sbi", p0_double_exp = [2000, 0.8, 0.1, 0.000001])
string += psg_vs_timing_without_bra_and_sbi
print(string)

# Write the results to a file
with open("two_term_exponential_distribution_results.tsv", "w") as f:
    f.write(string)