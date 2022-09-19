"""
Impose random perturbations to the kinetic models with physiologically relevant dynamics.
Prune them based on how fast they return to their original steady-state .
"""

from pytfa.io.json import load_json_model
from skimpy.io.yaml import load_yaml_model
from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, \
    load_concentrations, load_equilibrium_constants
from skimpy.analysis.ode import sample_initial_concentrations
from skimpy.core.solution import ODESolutionPopulation
from skimpy.core.parameters import ParameterValuePopulation, \
    load_parameter_population
from skimpy.utils.namespace import QSSA
from skimpy.utils.tabdict import TabDict
import os
import numpy as np
import pandas as pd
import multiprocessing as mp
import matplotlib.pyplot as plt

# Parameters
TIME = np.logspace(-9, np.log10(26.2), 1000)  # Iintegration time just above the doubling time of 26.1 hr
CONCENTRATION_SCALING = 1e6  # 1 mol to 1 mmol
FOLD_PERTURBATION = 2  # change in concentartions
NCPU = 10
N_SAMPLES = 50  # number of random perturbations imposed

# Paths
path_to_kmodel = '../../models/kinetic_models/kinetic_model.yaml'
path_to_tmodel = '../../models/tfa/tfa_model_no_interventions.json'
path_to_samples = './conc_samples/conc_samples_no_interventions.csv'
path_to_params = './output/pruned_parameters.hdf5'
path_to_output = './basins_of_attraction/ode_solutions_{}.csv'
path_for_pruned_params = './output/non_linear_pruned_parameters.hdf5'
path_to_plot = './basins_of_attraction/basins_of_attraction_results_sample_{}.png'
# Load pytfa model
tmodel = load_json_model(path_to_tmodel)
CPLEX = 'optlang-cplex'
tmodel.solver = CPLEX

tmodel.solver.configuration.tolerances.feasibility = 1e-9
tmodel.solver.configuration.tolerances.optimality = 1e-9
tmodel.solver.configuration.tolerances.integrality = 1e-9

# Load kinetic parameters
parameter_population = load_parameter_population(path_to_params)

# Load tfa samples
samples = pd.read_csv(path_to_samples, header=0, index_col=0)

# Compile ODE-Functions
kmodel = load_yaml_model(path_to_kmodel)
kmodel.prepare()
kmodel.compile_ode(sim_type=QSSA, ncpu=NCPU)
this_samples = list(parameter_population._index.keys())
samples_to_simulate = this_samples


def perturb_ode_solver(ix):
    print(ix)
    # Get reference concentrations
    tfa_id, _ = ix.split(',')
    tfa_id = int(tfa_id)
    sample = samples.loc[tfa_id]
    reference_concentrations = load_concentrations(sample, tmodel, kmodel,
                                                   concentration_scaling=CONCENTRATION_SCALING)
    # Add parameters
    kmodel.parameters = parameter_population[str(ix)]

    # Add a perturbation vector N_SAMPLES times to the reference concentrations
    # to create N_SAMPLES perturbed initial conditions
    perturbed_concentrations = sample_initial_concentrations(kmodel,
                                                             reference_concentrations,
                                                             lower_bound=1 / FOLD_PERTURBATION,
                                                             upper_bound=FOLD_PERTURBATION,
                                                             n_samples=N_SAMPLES)

    # Function to stop integration if the system diverges
    def rootfn(t, y, g, user_data):
        y0 = user_data['y0']
        n_max = user_data['max_norm']
        norm = np.sqrt(np.square((y - y0) / y0).sum())
        if norm >= n_max:
            g[0] = 0
        else:
            g[0] = 1

    user_data = {'y0': reference_concentrations[kmodel.variables].values,
                 'max_norm': 1e5}
    solutions = []

    # Loop through all perturbation vectors, and study their evolution in time
    for i, this_perturbations in perturbed_concentrations.iterrows():
        print(i)
        # Integrate perturbations
        kmodel.initial_conditions = TabDict([(k, v) for k, v in this_perturbations.iteritems()])
        this_sol_qssa = kmodel.solve_ode(TIME,
                                         solver_type='cvode',
                                         rtol=1e-9,
                                         atol=1e-9,
                                         max_steps=1e9,
                                         rootfn=rootfn,
                                         nr_rootfns=1,
                                         user_data=user_data
                                         )

        # Normalize the concentrations
        variables = this_sol_qssa.concentrations.columns
        this_sol_qssa.concentrations = (this_sol_qssa.concentrations - reference_concentrations[variables]) \
                                       / reference_concentrations[variables]
        solutions.append(this_sol_qssa)

    damped_time = -1  # calculate the damping effect after one cell cylce
    min_thresh = 0.05  # The dumping effect cutooff value
    save_solutions = 1
    results = []
    # If the model achieves 95% damping of the original perturbation then it is considered (non_linearly) stable
    for i in range(len(solutions)):
        final_norm = np.linalg.norm(solutions[i].concentrations.iloc[damped_time])
        original_norm = np.linalg.norm(solutions[i].concentrations.iloc[0])
        print(final_norm / original_norm)
        if final_norm > min_thresh * original_norm:
            results.append(0)
            save_solutions = 0
            # break
        else:
            results.append(1)

    if save_solutions == 1:
        mean_conc = solutions[0].concentrations
        for i in range(1, len(solutions)):
            mean_conc = mean_conc.add(solutions[i].concentrations)

        mean_conc = mean_conc / len(solutions)
        # plot
        plt.plot(solutions[0].time, mean_conc)
        plt.xlabel('time (h)')
        plt.ylabel('normalized conc')
        plt.savefig(path_to_plot.format(tfa_id))
        plt.close()

        solpop = ODESolutionPopulation(solutions)
        solpop.data.to_csv(path_to_output.format(ix))
        # Clear memory
        del solpop

    return results


# NOTE: maybe the above code must be run for one iteration before doing it in parallel
pool = mp.Pool(int(NCPU))
solutions = pool.starmap_async(perturb_ode_solver, [(ix,) for ix in samples_to_simulate])
pool.close()

solutions = solutions.get()
pool.terminate()

sturdy_params = []
sturdy_index = []
for i in range(0, len(solutions)):
    sol = solutions[i]
    if sum(sol) == N_SAMPLES:
        ind = list(parameter_population._index.keys())[i]
        sturdy_params.extend([parameter_population._data[i]])
        sturdy_index.extend([list(parameter_population._index.keys())[i]])

sturdy_parameters = ParameterValuePopulation(sturdy_params,
                                             kmodel=kmodel,
                                             index=sturdy_index)

sturdy_parameters.save(path_for_pruned_params)