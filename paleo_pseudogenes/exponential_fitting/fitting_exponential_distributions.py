import numpy as np
import lmfit
import pandas as pd

# Get dataframe with all the data
df = pd.read_csv("GeneralCharacteristicsSpecies.tsv", sep = "\t")
y_psg = df["WGM_pseudogenes"][1:] # remove Amborella from the list
y_ap = df["gene_anchorpairs_no_intron"][1:]/15
x_timing = df["most_recent_WGM(mya)"][1:]
x_Ks = df["most_recent_WGM(Ks_mean)"][1:].astype('float64')

def fit_curves(x, y, name, p0_exp = [2000, 0.02], p0_double_exp = [2000, 0.5, 0.01, 0.0001]):
    # Fit the data using lmfit: exponential
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
    print("R-squared (own calculation):", r_squared_exp_lmfit)

    ## Obtain parameters from the fitted curve
    a_exp = fit.best_values['a']
    b_exp = fit.best_values['b']

    chi_squared_exp = fit.chisqr
    reduced_chi_squared_exp = fit.redchi
    aic_exp = fit.aic
    bic_exp = fit.bic
    rsquared_exp = fit.rsquared
    rmse_exp = fit.residual.std()
    print("Chi-Squared:", chi_squared_exp)
    print("Reduced Chi-Squared:", reduced_chi_squared_exp)
    print("AIC:", aic_exp)
    print("BIC:", bic_exp)
    print("R-squared:", rsquared_exp)
    print("RMSE:", rmse_exp)

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
    print("R-squared (own calculation):", r_squared_double_exp_lmfit)

    ## Obtain parameters from the fitted curve
    A_double_exp = fit.best_values['A']
    q_double_exp = fit.best_values['q']
    a_double_exp = fit.best_values['a']
    c_double_exp = fit.best_values['c']

    chi_squared_d_exp = fit.chisqr
    reduced_chi_squared_d_exp = fit.redchi
    aic_d_exp = fit.aic
    bic_d_exp = fit.bic
    rsquared_d_exp = fit.rsquared
    rmse_d_exp = fit.residual.std()

    print("Chi-Squared:", chi_squared_d_exp)
    print("Reduced Chi-Squared:", reduced_chi_squared_d_exp)
    print("AIC:", aic_d_exp)
    print("BIC:", bic_d_exp)
    print("R-squared:", rsquared_d_exp)
    print("RMSE:", rmse_d_exp)

    result_string = name + "\t" + str(a_exp) + "\t" + str(b_exp) + "\t" +\
        str(chi_squared_exp) + "\t" + str(reduced_chi_squared_exp) + "\t" + str(aic_exp) + "\t" + str(bic_exp) + "\t" +\
        str(rsquared_exp) + "\t" + str(rmse_exp) + "\t" + str(A_double_exp) + "\t" + str(q_double_exp) + "\t" +\
        str(a_double_exp) + "\t" + str(c_double_exp) + "\t" +\
        str(chi_squared_d_exp) + "\t" + str(reduced_chi_squared_d_exp) + "\t" + str(aic_d_exp) + "\t" + str(bic_d_exp) + "\t" +\
        str(rsquared_d_exp) + "\t" + str(rmse_d_exp) + "\n"
    return result_string

string = "comparison\ta_exp\tb_exp\tchi_squared_exp\treduced_chi_squared_exp\taic_exp\tbic_exp\trsquared_exp\tRMSE_exp\tA_double_exp\tq_double_exp\ta_double_exp\tc_double_exp\tchi_squared_d_exp\treduced_chi_squared_d_exp\taic_d_exp\tbic_d_exp\trsquared_d_exp\tRMSE_d_exp\n"
psg_vs_Ks = fit_curves(x_Ks, y_psg, "psg_vs_Ks", p0_double_exp = [1000, 0.5, 0.1, 0.0001])
string += psg_vs_Ks
psg_vs_timing = fit_curves(x_timing, y_psg, "psg_vs_timing", p0_double_exp = [1000, 0.5, 0.1, 0.0001])
string += psg_vs_timing
APs_vs_Ks = fit_curves(x_Ks, y_ap, "APs_vs_Ks", p0_double_exp = [2000, 0.5, 0.01, 0.01])
string += APs_vs_Ks
APs_vs_timing = fit_curves(x_timing, y_ap, "APs_vs_timing", p0_double_exp = [2000, 0.5, 0.01, 0.001])
string += APs_vs_timing

# Write the results to a file
with open("curve_fitting_results.tsv", "w") as f:
    f.write(string)

# Too high (reduced) chi-squared values may be due to wrong estimate of SD because low data points (see: https://www.graphpad.com/support/faq/some-nonlinear-regression-programs-report-the-chi-square-of-a-fit-what-does-this-mean-why-doesnt-prism-report-the-chi-square-value/)
# AIC/BIC differences: https://psu-psychology.github.io/psy-597-SEM/09_model_comparison/model_comparison.html#model-evidence-k-l-distance (Burnham & Anderson, 2002)
# -> models in the 4–7 point range have considerably less support