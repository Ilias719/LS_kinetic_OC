"""
Based on a given population of kinetic parameters, modal analysis constructs pools of metabolites that are linked to the
system's eigenvalues. This code extracts the metabolites connected with slow eigenvalues.
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
import multiprocessing as mp
from skimpy.analysis.modal.modal_matrix import modal_matrix
from skimpy.core.parameters import ParameterValuePopulation, \
    load_parameter_population

NCPU = 28
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
TIMES_FASTER_THAN_MAX_EIG = 10  # Used for more strict pruning

# The range of steady-state samples used for the production of kinetic parameters
LOWER_IX = 0
UPPER_IX = 1000

# Number of stable sets of kinetic parameters per steady-state point
N_PARAMETER_SAMPLES_FOR_PARALLEL = 100

# Paths
path_to_kmodel = '../../models/kinetic_models/kinetic_model.yaml'
path_to_tmodel = '../../models/tfa/tfa_model_no_interventions.json'
path_to_samples = './conc_samples/conc_samples_no_interventions.csv'
path_to_params = './output/parameters_no_intervention_{}.hdf5'
path_to_results = './modal_analysis_results/modal_analysis_results.csv'

# Load the pytfa model
tmodel = load_json_model(path_to_tmodel)

CPLEX = 'optlang-cplex'
tmodel.solver = CPLEX

tmodel.solver.configuration.tolerances.feasibility = 1e-9
tmodel.solver.configuration.tolerances.optimality = 1e-9
tmodel.solver.configuration.tolerances.integrality = 1e-9
# tmodel.solver.problem.parameters.read.scale = -1
# tmodel.solver.problem.parameters.emphasis.numerical = 1
# tmodel.solver.configuration.presolve = True

# Load the kinetic model and prepare for compilinig of jacobian
kmodel = load_yaml_model(path_to_kmodel)
kmodel.prepare()
kmodel.compile_jacobian(ncpu=NCPU)
print('kinetic model jacobian --> compiled!')

samples = pd.read_csv(path_to_samples, header=0, index_col=0).iloc[LOWER_IX:UPPER_IX]

# Before doing it in parallel you have to initialize some stuff  in the kmodel, so we run it for one sample
stored_results = {name: [] for name, symbol in kmodel.variables.items()}
for i in range(LOWER_IX, LOWER_IX + 1):
    print('Steady state sample: ' + str(i))
    parameter_population = load_parameter_population(path_to_params.format(i))
    sample = samples.iloc[i, :]
    conc = load_concentrations(sample, tmodel, kmodel,
                               concentration_scaling=CONCENTRATION_SCALING)
    for j in range(0, 1):
        print(j)
        params = parameter_population[str(j)]
        params_with_strings = {str(name): value for name, value in params.items()}
        results = modal_matrix(kmodel, conc, params, params_with_strings)
        new_ind = [np.real(i) for i in results.index]
        results.index = new_ind
        sort_res = results.sort_index(ascending=False)
        slow_eig_ind = np.where(abs(sort_res.index) < TIMES_FASTER_THAN_MAX_EIG * abs(MAX_EIGENVALUES))
        for k in slow_eig_ind[0]:
            slow_met_ind = np.where(sort_res.iloc[k].abs() > 1e-2)
            temp = sort_res.iloc[k, slow_met_ind[0]]
            for ll in temp.keys():
                stored_results[ll].append(np.real(temp.loc[ll]))

# Make the modal analysis a function and then run in parallel
N_PARAMETER_SAMPLES = N_PARAMETER_SAMPLES_FOR_PARALLEL


def test_fun(sample_number):
    stored_results = {name: [] for name, symbol in kmodel.variables.items()}
    i = sample_number
    print('Steady state sample: ' + str(i))
    parameter_population = load_parameter_population(path_to_params.format(i))
    sample = samples.iloc[i, :]
    conc = load_concentrations(sample, tmodel, kmodel,
                               concentration_scaling=CONCENTRATION_SCALING)
    for j in range(0, N_PARAMETER_SAMPLES):
        print(j)
        params = parameter_population[str(j)]
        params_with_strings = {str(name): value for name, value in params.items()}
        results = modal_matrix(kmodel, conc, params, params_with_strings)
        new_ind = [np.real(i) for i in results.index]
        results.index = new_ind
        sort_res = results.sort_index(ascending=False)
        slow_eig_ind = np.where(abs(sort_res.index) < TIMES_FASTER_THAN_MAX_EIG * abs(MAX_EIGENVALUES))
        for k in slow_eig_ind[0]:
            slow_met_ind = np.where(sort_res.iloc[k].abs() > 1e-2)
            temp = sort_res.iloc[k, slow_met_ind[0]]
            for ll in temp.keys():
                stored_results[ll].append(np.real(temp.loc[ll]))
    return stored_results


# Run in parallel for faster results
pool = mp.Pool(int(25))
final_results = pool.starmap_async(test_fun, [(i,) for i in range(LOWER_IX, UPPER_IX)])
pool.close()
final_results = final_results.get()
pool.terminate()

# Save the results
stored_results = {name: [] for name, symbol in kmodel.variables.items()}
for i in range(LOWER_IX, UPPER_IX):
    for j in stored_results.keys():
        for ll in final_results[i][j]:
            stored_results[j].append(ll)

# Save the metabolite names, number of ocurrencies, number of reactions, average weight values
met_names = []
no_of_ocurr = []
total_param_samples = (UPPER_IX - LOWER_IX) * N_PARAMETER_SAMPLES
avg_value = []
pos_values = []
no_of_rxns = []
for i in stored_results.keys():
    if len(stored_results[i]) != 0:
        met_names.append(i)
        no_of_ocurr.append(len(stored_results[i]) / (460 * total_param_samples) * 100)
        avg_value.append(np.mean(abs(np.array(stored_results[i]))))
        pos_values.append(sum(np.array(stored_results[i]) > 0) / len(stored_results[i]))
        if i.startswith('_'):
            i = i[1:]
        no_of_rxns.append(len(list(tmodel.metabolites.get_by_id(i).reactions)))

stats = {'met_names': met_names, '% no of ocurr': no_of_ocurr, \
         'average abs weight': avg_value, '% of positive weights': pos_values, 'no of rxns': no_of_rxns}

df = pd.DataFrame(stats)
rounded = df.round({'% no of ocurr': 2, 'average abs weight': 4, '% of positive weights': 2})
rounded.sort_values('no of rxns', axis=0, inplace=True)
rounded.to_csv(path_to_results, index=False)