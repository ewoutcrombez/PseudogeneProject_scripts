# Figure plotting and fitting curves

These scripts were used to generate figures and fit the exponential functions to the time vs number data.

`exponential_fitting.py`:
To explore the evolution through time of WGM-derived pseudogenes and the number of gene-gene anchor pairs, the quantities were graphed against the timing of the most recent WGM, utilizing both million years ago (MYA) and Ks-value of the most recent WGM in each species (obtained from Vanneste et al., 2014). Subsequently, we attempted to fit a single- and double-exponential function to this dataset using the “curve_fit” function of the SciPy library (Virtanen et al., 2020) and using the lmfit library (Newville et al., 2016) in Python. The single-exponential function corresponds to the simple passive or random loss model whereby (pseudo)genes are lost randomly. The double-exponential function approximates the two-phase model described by Inoue et al. (2015). This model comprises an initial rapid phase, characterized by the simul-taneous loss of multiple genes, followed by a second phase with a constant rate of loss.
