"""
Code to sample populations of kinetic parameters based on the steady-state samples provdided.
The code has beem modified to run in parallel for faster sampling (the log messages have not been modified).
An additional pruning step saves all the kinetic parameters with physilogically relevant dynamics.
Finally, plots of distrubutions of the max-eigenvalues reported are created.
"""
from pytfa.io.json import load_json_model

from skimpy.io.yaml import load_yaml_model
from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, \
    load_concentrations, load_equilibrium_constants
from skimpy.core.parameters import ParameterValuePopulation, load_parameter_population
from skimpy.sampling.simple_parameter_sampler import SimpleParameterSampler
from skimpy.utils.general import sanitize_cobra_vars
from scipy.sparse import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp

# The range of steady-state samples that will be used for the production of kinetic parameters
LOWER_IX = 0
UPPER_IX = 1000

NCPU = 25
# Number of stable sets of kinetic parameters per steady-state point
N_PARAMETER_SAMPLES_FOR_PARALLEL = 100

# Scaling parameters
CONCENTRATION_SCALING = 1e6  # 1 mol to 1 mmol
TIME_SCALING = 1  # 1hour

# Parameters of human cancer cells
DENSITY = 1200  # g/L
GDW_GWW_RATIO = 0.25  # Assumes 75% Water

flux_scaling_factor = 1e-3 * (GDW_GWW_RATIO * DENSITY) * CONCENTRATION_SCALING / TIME_SCALING

# Doubling time 26.1 hours from the biomass upper bound (ln(2)/0.026)
MAX_EIGENVALUES = -3 * 1 / 26.1  # 3 times faster than the doubling time of the cell (26.1 hours)
MIN_EIGENVALUES = -1 / 2.7e-10  # 1/1mus = 1/1e-6/60/60 1/h

# Paths
path_to_kmodel = '../../models/kinetic_models/kinetic_model.yaml'
path_to_tmodel = '../../models/tfa/tfa_model_no_interventions.json'
path_to_samples = './conc_samples/conc_samples_no_interventions.csv'
path_to_lambda_values = './max_eigenvalues_no_intervention_{}_{}.csv'
path_for_output = './output/parameters_no_intervention_{}.hdf5'
path_for_params = './output/pruned_parameters_.hdf5'
path_to_plot = './graphs/max_eignevalues_no_intervention.png'
# Load pytfa model
tmodel = load_json_model(path_to_tmodel)

CPLEX = 'optlang-cplex'
tmodel.solver = CPLEX

tmodel.solver.configuration.tolerances.feasibility = 1e-9
tmodel.solver.configuration.tolerances.optimality = 1e-9
tmodel.solver.configuration.tolerances.integrality = 1e-9
# tmodel.solver.problem.parameters.read.scale = -1
# tmodel.solver.problem.parameters.emphasis.numerical = 1
# tmodel.solver.configuration.presolve = True

# Load steady-state samples
samples = pd.read_csv(path_to_samples, header=0, index_col=0).iloc[LOWER_IX:UPPER_IX]

# Load the kinetic model and prepare for compilinig of jacobian
kmodel = load_yaml_model(path_to_kmodel)
kmodel.prepare()
kmodel.compile_jacobian(ncpu=NCPU)
print('kinetic model jacobian --> compiled!')


# Define the kinetic parameter sampling function
def kin_sample(i, sample, tmodel_=tmodel, kmodel_=kmodel):
    print('Producing kinetic models from sample:', i)
    sampler_params = SimpleParameterSampler.Parameters(n_samples=N_PARAMETER_SAMPLES)
    sampler = SimpleParameterSampler(sampler_params)

    # Load fluxes and concentrations
    fluxes = load_fluxes(sample, tmodel_, kmodel_,
                         density=DENSITY,
                         ratio_gdw_gww=GDW_GWW_RATIO,
                         concentration_scaling=CONCENTRATION_SCALING,
                         time_scaling=TIME_SCALING)

    concentrations = load_concentrations(sample, tmodel_, kmodel_,
                                         concentration_scaling=CONCENTRATION_SCALING)

    # Fetch equilibrium constants
    load_equilibrium_constants(sample, tmodel_, kmodel_,
                               concentration_scaling=CONCENTRATION_SCALING,
                               in_place=True)

    params, lamda_max, lamda_min = sampler.sample(kmodel_,
                                                  fluxes,
                                                  concentrations,
                                                  only_stable=True,
                                                  min_max_eigenvalues=True)
    params_population = ParameterValuePopulation(params, kmodel_)
    params_population.save(path_for_output.format(i))
    return [lamda_max, lamda_min]


# Due to some class initializations, it is mandatory to do a first sampling run not in parallel
i = LOWER_IX
N_PARAMETER_SAMPLES = 10
sample = samples.iloc[0, :]
lambda_max_min = kin_sample(i, sample)

# Run the process in parallel
N_PARAMETER_SAMPLES = N_PARAMETER_SAMPLES_FOR_PARALLEL
pool = mp.Pool(int(NCPU))
lambda_max_min = pool.starmap_async(kin_sample, [(i, sample) for i, sample in samples.iterrows()])
pool.close()
lambda_max_min = lambda_max_min.get()
lambda_min = {}
lambda_max = {}
pool.terminate()

# Curate the results
for i in range(0, len(lambda_max_min)):
    lambda_max['sample ' + str(i)] = lambda_max_min[i][0]
    lambda_min['sample ' + str(i)] = lambda_max_min[i][1]

lambda_df = pd.DataFrame.from_dict(lambda_max, orient='index')

lambda_df.to_csv(path_to_lambda_values.format(LOWER_IX, UPPER_IX))

# Prune parameters based on the eigenvalue cutoffs
is_selected = (lambda_df < MAX_EIGENVALUES)  # & (lambda_min_all > MIN_EIGENVALUES )
is_selected.columns = range(lambda_df.shape[1])
is_selected.index = samples.index
fast_parameters = []
fast_index = []

for i, row in is_selected.iterrows():
    if any(row):
        print(i)
        fast_models = np.where(np.array(row))[0]

        # Load the respective solutions
        parameter_population = load_parameter_population(path_for_output.format(i))

        # Construct the Parameter populations for each sample
        fast_parameters.extend([parameter_population._data[k] for k in fast_models])
        fast_index.extend(["{},{}".format(i, k) for k in fast_models])

print("{} models found for {} < lam < {}".format(len(fast_index),
                                                 MIN_EIGENVALUES,
                                                 MAX_EIGENVALUES))
# Save the fast models
fast_parameters = ParameterValuePopulation(fast_parameters,
                                           kmodel=kmodel,
                                           index=fast_index)

fast_parameters.save(path_for_params)

lambda_max = pd.read_csv(path_to_lambda_values.format(LOWER_IX, UPPER_IX), header=0, index_col=0)
ar = np.array(lambda_max).reshape((UPPER_IX - LOWER_IX) * N_PARAMETER_SAMPLES, 1)

fig = plt.figure(figsize=(10, 7))
plt.boxplot(-ar, whis=1e6)  # need to have positive values in order to assign log scale
plt.yscale('log')
# plt.plot(1, -MAX_EIGENVALUES, marker='o', markersize=10)
plt.axhline(y=-MAX_EIGENVALUES, color='r', linestyle='-')
plt.xticks(np.ceil([N_PARAMETER_SAMPLES / 2]), ['Max eigenvalues no intervention'])
plt.ylim([1e-9, 1])
plt.yticks(np.logspace(-9, 0, 10))
plt.savefig(path_to_plot)
plt.close()

# Find the k biggest max eigenvalues
# sorted_ar = np.sort(ar.T)
# k = 15
# fast_sample_id = []
# for i in range(0, k):
#     pos = np.where(lambda_max == sorted_ar[0][i])
#     fast_sample_id.append(pos[0][0])
