"""
Code that estimates possible drug targtes based on Metabolic Control Analysis.
"""
from pytfa.io.json import load_json_model
from skimpy.io.yaml import load_yaml_model
from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, \
    load_concentrations, load_equilibrium_constants
from skimpy.analysis.ode import sample_initial_concentrations
from skimpy.core.solution import ODESolutionPopulation
from skimpy.core.parameters import ParameterValuePopulation, \
    load_parameter_population
from skimpy.analysis.mca.utils import get_dep_indep_vars_from_basis
from skimpy.utils.general import get_stoichiometry
from skimpy.utils.tensor import Tensor

from skimpy.utils.namespace import QSSA
from skimpy.utils.tabdict import TabDict

import numpy as np
import pandas as pd

import os.path

# Parameters
NCPU = 28

# kmodel.parameters.items use the following name. However, the actual parameter is vmax_forward_ so watch out on
# transport_reactions!
PARAMETER_FOR_MCA = 'vmax_forward'

# Parameters
CONCENTRATION_SCALING = 1e6  # 1 mol to 1 mumol
TIME_SCALING = 1  # 1hr to 1min
DENSITY = 1200  # g/L
GDW_GWW_RATIO = 0.25  # Assumes 75% Water

# Paths
path_to_kmodel = '../../models/kinetic_models/kinetic_model_18_7.yaml'
path_to_tmodel = '../../models/tfa_model_ATP_modified_modal_analysis_original_code_1e-3_flux_fast_linear_all.json'
path_to_tfa_samples = '../kinetics/conc_samples/conc_samples_ATP_modified_10_fastest.csv'
path_to_params = '../kinetics/output/non_linear_pruned_parameters.hdf5'
path_to_output = "./output/basins/ode_solutions_{}.csv"

# load models
tmodel = load_json_model(path_to_tmodel)
kmodel = load_yaml_model(path_to_kmodel)

parameter_population = load_parameter_population(path_to_params)

# Load tfa samples
tfa_samples = pd.read_csv(path_to_tfa_samples, header=0, index_col=0)

# Prepare the kinetic model
kmodel.prepare()

# Compile with parameter elasticities
parameter_list = TabDict([(k, p.symbol) for k, p in kmodel.parameters.items()
                          if p.name.startswith(PARAMETER_FOR_MCA)])

kmodel.compile_mca(sim_type=QSSA, ncpu=NCPU, parameter_list=parameter_list)

this_samples = list(parameter_population._index.keys())
samples_to_simulate = this_samples

flux_control_data = []

# Take only the samples that passed the nonlinear response test
nonlinear_samples = samples_to_simulate

for ix in nonlinear_samples:
    print(ix)
    # Add the kinetic parameters of the current set to the kinetic model
    kmodel.parameters = parameter_population[ix]

    # Get the TFA sample corresponding to the current kinetic parameter set
    tfa_id, _ = ix.split(',')
    tfa_id = int(tfa_id)
    tfa_sample = tfa_samples.loc[tfa_id]

    # Load steady state fluxes and concentrations from the TFA sample into the kinetic model
    fluxes = load_fluxes(tfa_sample, tmodel, kmodel,
                         density=DENSITY,
                         ratio_gdw_gww=GDW_GWW_RATIO,
                         concentration_scaling=CONCENTRATION_SCALING,
                         time_scaling=TIME_SCALING)

    concentrations = load_concentrations(tfa_sample, tmodel, kmodel,
                                         concentration_scaling=CONCENTRATION_SCALING)

    flux_control_coeff = kmodel.flux_control_fun(fluxes, concentrations, [parameter_population[ix]])
    flux_control_data.append(flux_control_coeff._data)

# Make tensor for the population level
fcc_data = np.concatenate(flux_control_data, axis=2)

flux_index = pd.Index(kmodel.reactions.keys(), name="flux")
parameter_index = pd.Index(kmodel.flux_control_fun.parameter_elasticity_function.respective_variables, name="parameter")
sample_index = pd.Index(nonlinear_samples, name="sample")

# Get the mean flux control coefficients for this sample
flux_control_coeff = Tensor(fcc_data, [flux_index, parameter_index, sample_index])
mean_fcc = flux_control_coeff.mean('sample')

# Remove transports

transport_reactions =  [PARAMETER_FOR_MCA+'_'+r.id for r in tmodel.reactions
                        if (len(r.compartments) > 1 )
                        and not ('i' in r.compartments)
                        and not r.id.startswith('LMPD_')
                        and r.id in kmodel.reactions]

mean_fcc = mean_fcc.drop(transport_reactions, axis=1)


# TOP20 Control coefficents for phenylalanine overproduction
sorted_index = mean_fcc.loc['biomass'].abs().sort_values().index[-20:]
fcc_biomass_slice = flux_control_coeff.slice_by('flux', 'biomass')

# Bar plot
import matplotlib.pyplot as plt
f = plt.figure(figsize=(4,8))
ax = f.add_subplot(111)

plt.barh(range(len(sorted_index)),  fcc_biomass_slice.loc[sorted_index, :].mean(axis=1),)
plt.xlabel(r'Control for biomass production')

''' 
Make sure you have the correct tick labels

plt.yticks( _____ ) 
'''
plt.yticks(range(0, 20), [(sorted_index[i])[13::] for i in range(len(sorted_index))])
plt.tight_layout()
plt.savefig('./bars_control_coeff.png')

plt.close()



