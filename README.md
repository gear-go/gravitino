# gravitino
Exploration of the SUSY parameter space with MultiNest to look for points compatible with the cosmic ray lepton measurements

Gravitino_multinest.ipynb: An example of using MultiNest to fit a phenomenological model. It shows how to upload the data, create priors, define likelihood, run MultiNest and plot results.

The different codes to fit gravitino models to data are:
fit_all_leptons_new_data_multinest_AMS02Sum_gravitino_shade.py

fit_all_leptons_new_data_multinest_CALET_shade.py

fit_all_leptons_new_data_multinest_dampe_gravitino_shade.py

fit_all_leptons_new_data_multinest_noSum_gravitino_shade.py

After running one of the files you will produce a series of valid points. To marginalize to find confidence regions we use the Johannes' code "multinest_marginals.py". It runs as follows using the prefix used in the fitting code:

 python multinest_marginals.py fit_multines_C_Gravitino_DAMPE_shade_4_
 
 The output of that code represents the favourable regions of the parameter space that can fit the data.
