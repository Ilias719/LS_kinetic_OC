"""
Based on the provided mat file, the code extracts the FDP and the bounds of all the variables (LC, DG, ...).
It then imposes basal fluxes, min thermodynamic displacement and forces diffusion transport reactions to operate close
to equilibrium.
Lastly, the code does variability analysis and runs an initial sampling of steady-state sample points.
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
from pytfa.io.json import save_json_model
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

# Path to data and model
path_to_mat_tmodel = './../../models/tfa/model_remi.mat'
path_to_thermodb = './../../models/tfa/recon3.thermodb'
path_to_tmodel =  '../../models/tfa/tfa_model_no_interventions.json'
path_to_samples = './conc_samples/conc_samples_no_interventions.csv'

# Load cobra model from mat file
cmodel = import_matlab_model(path_to_mat_tmodel)

# Convert to a thermodynamics model
thermo_data = load_thermoDB(path_to_thermodb)
tmodel = pytfa.ThermoModel(thermo_data, cmodel)

# Setup solver
CPLEX = 'optlang-cplex'
tmodel.solver = CPLEX

tmodel.solver.configuration.tolerances.feasibility = 1e-9
tmodel.solver.configuration.tolerances.optimality = 1e-9
tmodel.solver.configuration.tolerances.integrality = 1e-9
tmodel.solver.problem.parameters.emphasis.numerical = 1
# tmodel.solver.configuration.presolve = True
# tmodel.solver.problem.parameters.read.scale = -1

# TFA conversion
tmodel.prepare()
tmodel.convert(add_displacement=True)
solution = tmodel.optimize()

# Import variable bounds from mat file
mat_data = loadmat(path_to_mat_tmodel)
mat_model = mat_data['model_fdp'][0, 0]
varname_index = {e[0][0]: i for i, e in enumerate(mat_model['varNames'])}

# Extract DGo values from the mat file
for dgo in tmodel.delta_gstd:
    var_name = dgo.name
    i = varname_index[var_name]

    lb = mat_model['var_lb'][i][0]
    ub = mat_model['var_ub'][i][0]
    if dgo.variable.ub < lb:
        dgo.variable.ub = ub
        dgo.variable.lb = lb
    else:
        dgo.variable.lb = lb
        dgo.variable.ub = ub

# Extract LogConcentration values from the mat file
for lc in tmodel.log_concentration:
    var_name = lc.name
    if var_name in varname_index:
        i = varname_index[var_name]
    else:
        continue

    lb = mat_model['var_lb'][i][0]
    ub = mat_model['var_ub'][i][0]
    if lc.variable.ub < lb:
        lc.variable.ub = ub
        lc.variable.lb = lb
    else:
        lc.variable.lb = lb
        lc.variable.ub = ub

# Extract net-flux values from the mat file
for rxn in tmodel.reactions:
    forward_variable = rxn.forward_variable
    reverse_variable = rxn.reverse_variable

    i = varname_index['NF_{}'.format(rxn.id)]

    lb = mat_model['var_lb'][i][0]
    ub = mat_model['var_ub'][i][0]
    rxn.bounds = (lb, ub)

# Remove blocked reactions
blocked = ['r0122', 'PPAm']
tmodel.remove_reactions(blocked)

# Impose basal flux constraints
MIN_FLUX = 1e-6
add_min_flux_requirements(tmodel, MIN_FLUX, inplace=True)

solution = tmodel.optimize()
tmodel.reactions.biomass.lower_bound = 0.8 * solution.objective_value

# Add all missing DGs
tmodel = add_undefined_delta_g(tmodel,
                               solution,
                               delta_g_std=0.0,
                               delta_g_std_err=10.0,
                               add_displacement=True,
                               inplace=True)  # NOTE: exchange reactions do not have themodynamic properties

# Force a minimal thermodynamic displacement
tmodel_fdp = tmodel.copy()
tva_fluxes = variability_analysis(tmodel_fdp, kind='reactions')

min_log_displacement = 1e-3
tmodel_fdp = add_min_log_displacement(tmodel_fdp,
                                      min_log_displacement,
                                      tva_fluxes=tva_fluxes,
                                      inplace=False)

# Find the reactions that are diffusion transporters
transport_reactions = []
for i in range(0, len(mat_model['rxns'])):
    if mat_model['isTrans'][i] == 1:
        try:
            if len(tmodel_fdp.reactions.get_by_id(mat_model['rxns'][i][0][0]).metabolites) == 2:  # only diffusion
                transport_reactions.append(tmodel.reactions.get_by_id(mat_model['rxns'][i][0][0]).id)
        except:  # not needed
            continue

thermo_disp_varnames = [tmodel_fdp.thermo_displacement.get_by_id(r_id).name for r_id in transport_reactions]
tva_disp = [element_variability_analysis(tmodel_fdp, i) for i in thermo_disp_varnames]
tva_disp = pd.DataFrame(tva_disp, index=thermo_disp_varnames)

MAX_LOG_TRANSPORT_DISPLACEMENT = np.log(0.95)
for r_id in transport_reactions:
    variable = tmodel_fdp.thermo_displacement.get_by_id(r_id).variable
    min_disp, max_disp = tva_disp.loc[variable.name, :]

    ub_old = variable.ub
    lb_old = variable.lb

    print(max_disp, min_disp)

    if -MAX_LOG_TRANSPORT_DISPLACEMENT > min_disp > 0:
        print("UB: {}, NEW {} MIN {} MAX {}".format(variable.ub, -MAX_LOG_TRANSPORT_DISPLACEMENT, min_disp, max_disp))
        variable.ub = -MAX_LOG_TRANSPORT_DISPLACEMENT

    elif MAX_LOG_TRANSPORT_DISPLACEMENT < max_disp < 0:
        print("UB: {}, NEW {} MIN {} MAX {}".format(variable.lb, MAX_LOG_TRANSPORT_DISPLACEMENT, min_disp, max_disp))
        variable.lb = MAX_LOG_TRANSPORT_DISPLACEMENT
    else:
        print('Reaction needs to be displaced from eqiulibrium {}'.format(r_id))
        continue
    try:
        tmodel_fdp.optimize()
    except:
        print('Reaction needs to be displaced from eqiulibrium {}'.format(r_id))
        variable.ub = ub_old
        variable.lb = lb_old

# Do variability analysis and impose the new bounds
tva_fluxes = variability_analysis(tmodel_fdp, kind='reactions')
thermo_vars = [DeltaG, DeltaGstd, ThermoDisplacement, LogConcentration]
tva_thermo = variability_analysis(tmodel_fdp, kind=thermo_vars)
tight_model = apply_reaction_variability(tmodel_fdp, tva_fluxes)
tight_model = apply_generic_variability(tight_model, tva_thermo)
tight_model.objective = tight_model.reactions.biomass
sol = tight_model.optimize()

# Before moving on you might want to check for possible slow turnover ratios

# Remove integer variables in preparation for sampling
continuous_model = strip_from_integer_variables(tight_model)

# Do a first sampling run
N_SAMPLES = 1000
samples = sample(continuous_model, N_SAMPLES, method='achr', thinning=500)
samples.to_csv(path_to_samples)

# Save curated version of the model
save_json_model(tight_model, path_to_tmodel)


