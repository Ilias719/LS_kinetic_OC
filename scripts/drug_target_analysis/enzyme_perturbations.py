"""
This code is used to apply enzyme perturbations to the system and solve the ODEs in order to study the effects on
biomass production. Used both for drug target analysis and application of MCA results.
"""

from pytfa.io.json import load_json_model
from skimpy.io.yaml import load_yaml_model
from skimpy.analysis.oracle.load_pytfa_solution import load_concentrations
from skimpy.analysis.ode import sample_initial_concentrations
from skimpy.core.solution import ODESolutionPopulation
from skimpy.core.parameters import ParameterValuePopulation, \
    load_parameter_population
from skimpy.analysis.ode.utils import make_flux_fun
from skimpy.utils.namespace import QSSA
from skimpy.utils.tabdict import TabDict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os.path

NCPU = 28
TIME = np.logspace(-9, np.log10(52.2), 1000)  # Integration time: two times the doubling time of the cell (26.1 h)

# Scaling parameters
CONCENTRATION_SCALING = 1e6  # 1 mol to 1 mumol
TIME_SCALING = 1  # 1 hour to 1 min
DENSITY = 1200  # g/L
GDW_GWW_RATIO = 0.25  # Assumes 75% Water
flux_scaling_factor = 1e-3 * (GDW_GWW_RATIO * DENSITY) * CONCENTRATION_SCALING / TIME_SCALING


# Enzyme to perturb (1 or more) and the desired fold change for each of them
ENZYMES_TO_PERTURB = ['FOLOAT1tc', 'FOLOAT2tc', 'FOLABCCte', 'FOLOATPtc', 'r0963', 'FOLt2', 'FOLTle', '5MTHFt']
CHANGE =             [0, 0, 0, 0, 0, 0, 0, 0]

# Paths
path_to_kmodel = '../../models/kinetic_models/kinetic_model.yaml'
path_to_tmodel = '../../models/tfa/tfa_model_no_interventions.json'
path_to_samples = './conc_samples/conc_samples_no_interventions.csv'
path_to_params = '../kinetics/output/non_linear_pruned_parameters.hdf5'
path_to_output = './enzyme_perturbations_results/ode_solutions_folate_receptors_0%_{}.csv'
path_to_plots = './enzyme_perturbations_results/biomass_flux_folate_receptors_0%.png'#.format(ENZYMES_TO_PERTURB[0], CHANGE[0])

# load the pytfa model
tmodel = load_json_model(path_to_tmodel)

parameter_population = load_parameter_population(path_to_params)

# Load tfa samples
samples = pd.read_csv(path_to_samples, header=0, index_col=0)

# Load the kinetic model, prepare and compile ODEs
kmodel = load_yaml_model(path_to_kmodel)
kmodel.prepare(mca=False)
kmodel.compile_ode(sim_type=QSSA, ncpu=NCPU)

# Create enzmye perturbations list
ENZYME_PERTURBATIONS = [(str(kmodel.reactions[r_id].parameters.vmax_forward.symbol), change)
                        for r_id, change in zip(ENZYMES_TO_PERTURB, CHANGE)]

# Flux expressions based on odes
calc_flux = make_flux_fun(kmodel, QSSA)

this_samples = list(parameter_population._index.keys())
samples_to_simulate = this_samples
solutions = []

# For every set of kinetic parameters, apply the enzyme perturbation and integrate the system
for ix in samples_to_simulate:
    print(ix)

    # Add parameters
    kmodel.parameters = parameter_population[ix]

    # Get reference concentrations
    tfa_id, _ = ix.split(',')
    tfa_id = int(tfa_id)
    sample = samples.loc[tfa_id]
    concentrations = load_concentrations(sample, tmodel, kmodel,
                                                   concentration_scaling=CONCENTRATION_SCALING)

    # Here decide whether you want to apply all perturbations together or examine them one by one
    # Currently, this is set up to do all together
    for parameter, change in ENZYME_PERTURBATIONS:
        # Integrate perturbations
        print(min(concentrations), max(concentrations))
        kmodel.initial_conditions = TabDict([(k, v) for k, v in concentrations.iteritems()])

        # Perturb parameter
        kmodel.parameters[parameter].value = kmodel.parameters[parameter].value*change

    this_sol_qssa = kmodel.solve_ode(TIME,
                                        solver_type='cvode',
                                        rtol=1e-9,
                                        atol=1e-9,
                                        max_steps=1e9)

    solutions.append(this_sol_qssa)

# Save the solutions
solpop = ODESolutionPopulation(solutions, [str(x) for x in ENZYME_PERTURBATIONS])
solpop.data.to_csv(path_to_output.format(ix))

# Clear memory
del solpop

# Calculate biomass flux values
biomass_flux = {}
for i in range(0, len(solutions)):
    print(i)
    tmp_bio = []
    params = {k: v for k, v in parameter_population[samples_to_simulate[i]].items()}
    for kk in range(0, solutions[i].concentrations.shape[0]):
        conc = solutions[i].concentrations.iloc[kk]
        fluxes = pd.Series(calc_flux(conc, parameters=params))
        tmp_bio.append(fluxes.loc['biomass'] / flux_scaling_factor)
    biomass_flux[samples_to_simulate[i]] = tmp_bio


# Plot all the results
for i in biomass_flux.keys():
    plt.plot(solutions[0].time, biomass_flux[i])
plt.tight_layout()
plt.savefig(path_to_plots)
plt.close()

