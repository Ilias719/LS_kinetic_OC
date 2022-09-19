"""
Code that samples steady-state sample points
"""
import numpy as np
import pandas as pd
from scipy.io import loadmat
import pytfa
from pytfa.io import import_matlab_model, load_thermoDB
from pytfa.analysis import variability_analysis, \
    apply_reaction_variability, \
    apply_generic_variability, \
    apply_directionality
from pytfa.optim.variables import DeltaG, DeltaGstd, ThermoDisplacement, LogConcentration
from pytfa.io.json import save_json_model, load_json_model
from pytfa.optim.relaxation import relax_dgo
from skimpy.analysis.oracle import *
from skimpy.core.compartments import Compartment
from skimpy.utils.general import sanitize_cobra_vars
from sys import argv
from skimpy.io.generate_from_pytfa import FromPyTFA
from skimpy.io.yaml import export_to_yaml
from pytfa.optim import strip_from_integer_variables
from pytfa.analysis import sample
from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, \
    load_concentrations, load_equilibrium_constants
from skimpy.core.parameters import ParameterValuePopulation, load_parameter_population
from skimpy.sampling.simple_parameter_sampler import SimpleParameterSampler
from skimpy.utils.general import sanitize_cobra_vars
from scipy.sparse import *
from element_variability_analysis import  element_variability_analysis

# Paths
path_to_tmodel =  '../../models/tfa/tfa_model_no_interventions.json'
path_to_samples = './conc_samples/conc_samples_no_interventions.csv'


tmodel = load_json_model(path_to_tmodel)

# Setup solver
tmodel.solver.configuration.tolerances.feasibility = 1e-9
tmodel.solver.configuration.tolerances.optimality = 1e-9
tmodel.solver.configuration.tolerances.integrality = 1e-9
# tmodel.solver.problem.parameters.read.scale = -1
tmodel.solver.problem.parameters.emphasis.numerical = 1
# tmodel.solver.configuration.presolve = True

# Parameters
MIN_MAX = False
N_SAMPLES = 1000
thinning_factor = 500

# If there has been changes to the thermo model do a variability analysis to impose new bounds
if MIN_MAX:
    # Do variability analysis and impose the new bounds
    tva_fluxes = variability_analysis(tmodel, kind='reactions')
    thermo_vars = [DeltaG, DeltaGstd, ThermoDisplacement, LogConcentration]
    tva_thermo = variability_analysis(tmodel, kind=thermo_vars)
    tight_model = apply_reaction_variability(tmodel, tva_fluxes)
    tight_model = apply_generic_variability(tight_model, tva_thermo)
    tight_model.objective = tight_model.reactions.biomass
    sol = tight_model.optimize()

    # Remove integer variables in preparation for sampling
    continuous_model = strip_from_integer_variables(tight_model)

    # Sampling
    samples = sample(continuous_model, N_SAMPLES, method='achr', thinning=thinning_factor)

    # Save results
    samples.to_csv(path_to_samples)
    save_json_model(tight_model, path_to_tmodel[:-5]+'tight_json')
else:
    # Remove integer variables in preparation for sampling
    continuous_model = strip_from_integer_variables(tmodel)

    # Sampling
    samples = sample(continuous_model, N_SAMPLES, method='achr', thinning=thinning_factor)

    # Save results
    samples.to_csv(path_to_samples)
