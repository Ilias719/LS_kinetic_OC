"""
THIS CODE IS NOT NECESSARY

Code that calculates turnover ratios (sum(consuming fluxes)/(metabolite concentration)) based on a given set of
steady-state samples.

"""

from skimpy.analysis.oracle import *  # Important to avoid MinFluxVariable serialization error
from pytfa.io.json import load_json_model
from skimpy.io.yaml import load_yaml_model
from skimpy.utils.general import sanitize_cobra_vars
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sys import argv
from pytfa.analysis.variability import _variability_analysis_element

path_to_tmodel = '../../models/tfa/tfa_model_no_interventions.json'
path_to_samples = './conc_samples/conc_samples_no_interventions.csv'
path_to_kmodel = '../../models/kinetic_models/kinetic_model.yaml'

# Load pytfa model
tmodel = load_json_model(path_to_tmodel)
tmodel.solver.configuration.tolerances.feasibility = 1e-9
tmodel.solver.configuration.tolerances.optimality = 1e-9
tmodel.solver.configuration.tolerances.integrality = 1e-9

# Load the kinetic model
kmodel = load_yaml_model(path_to_kmodel)

tfa_samples = pd.read_csv(path_to_samples, index_col=0, header=0)

variable_mets = list(kmodel.reactants.keys())

# Scaling parameters
CONCENTRATION_SCALING = 1e6  # 1 mol to 1 mmol
TIME_SCALING = 1  # 1hour

# Parameters of human cancer cells
DENSITY = 1200  # g/L
GDW_GWW_RATIO = 0.25  # Assumes 75% Water

flux_scaling_factor = 1e-3 * (GDW_GWW_RATIO * DENSITY) * CONCENTRATION_SCALING / TIME_SCALING
eigenvalue_cutoff = 3 / 26.1  # 3 times faster than the doubling time of the cell (26.1 hours)

# For each steady-state sample, calculate the turnover ratios of the metabolites
min_turnover_metabolite = []
list_of_all_mets_slow_turnover = []
turnover = []
for i, solution in tfa_samples.iterrows():
    print('Sample '+str(i))
    this_solutions_turnover = []
    list_of_vars = []
    for this_lc in tmodel.log_concentration:

        if sanitize_cobra_vars(this_lc.id) in variable_mets:
            list_of_vars.append(this_lc.id)
            LC = this_lc.variable.name
            concentration = np.exp(solution[LC]) * CONCENTRATION_SCALING
            # For each metabolite calculate the metabolite concentration over the sum of all negative fluxes
            this_met = tmodel.metabolites.get_by_id(this_lc.id)
            rxns = this_met.reactions
            consuming_flux = 0
            for this_rxn in rxns:
                this_flux = (solution[this_rxn.forward_variable.name] - solution[this_rxn.reverse_variable.name]) \
                            * flux_scaling_factor
                # Test with outgoing flux
                if this_flux * this_rxn.metabolites[this_met] < 0:
                    consuming_flux += this_flux * this_rxn.metabolites[this_met]
            this_turnover = -consuming_flux / concentration
            this_solutions_turnover.append(this_turnover)

    this_solution_min_turnover_metabolite_index = \
        np.where(np.asarray(this_solutions_turnover) == np.min(this_solutions_turnover))[0][0]
    min_turnover_metabolite.append(
        [list_of_vars[this_solution_min_turnover_metabolite_index], np.min(this_solutions_turnover)])
    turnover.append(this_solutions_turnover)

    # Print those mets whose turnover is < eigenvalue_cutoff
    slow_turnover_mets = {list_of_vars[i]: this_solutions_turnover[i] for i in
                          list(np.where(np.asarray(this_solutions_turnover) < eigenvalue_cutoff)[0])}
    flag = 0
    flag = [k for k in slow_turnover_mets.keys() if len(tmodel.metabolites.get_by_id(k).reactions) == 2]
    if flag != []:
        print('This sample has a slow TC for a linear metabolite!')
        print(flag)
    list_of_all_mets_slow_turnover.append(slow_turnover_mets)

df_mets_slow = pd.DataFrame(list_of_all_mets_slow_turnover)
mean_turnovers = df_mets_slow.mean()

# Look for those mets in a linear pathway
dict_slow_turnovers = {}
for this_met in mean_turnovers.index:
    n_rxns = len(tmodel.metabolites.get_by_id(this_met).reactions)
    n_occurences = df_mets_slow.count().loc[this_met]
    exp_lb = np.exp(tmodel.log_concentration.get_by_id(this_met).variable.lb)
    exp_ub = np.exp(tmodel.log_concentration.get_by_id(this_met).variable.ub)
    dict_temp = {'n_reactions': n_rxns,
                 'n_occurences': n_occurences,
                 'mean_turnover': mean_turnovers.loc[this_met],
                 'lower_bound': exp_lb,
                 'upper_bound': exp_ub
                 }
    dict_slow_turnovers[this_met] = dict_temp

df_slow_turnovers = pd.DataFrame(dict_slow_turnovers)
df_slow_turnovers.T.sort_values(by='n_occurences', ascending=False).to_csv('slow_metabolites_information.csv')

# Add constraints to the problematic mets
lc_ub = np.log(1e-4)
failed_count = 0
failed_mets = []
for met_name, values in df_slow_turnovers.T.iterrows():
    if values.n_occurences > 0.1 * len(tfa_samples):
        # check if we can enforce concnentration constraint based on the flux lower bound
        list_reac = list(tmodel.metabolites.get_by_id(met_name).reactions)
        if np.abs(list_reac[0].lower_bound) > np.abs(list_reac[0].upper_bound):
            flux_value = list_reac[0].upper_bound
        else:
            flux_value = list_reac[0].lower_bound
        flux_lb = flux_value * list_reac[0].metabolites[tmodel.metabolites.get_by_id(met_name)] * flux_scaling_factor
        conc_ub = np.log(np.abs(flux_lb) / (5 * eigenvalue_cutoff) / CONCENTRATION_SCALING)
        print(met_name + ' ' + str(np.exp(conc_ub)))
        original_ub = tmodel.log_concentration.get_by_id(met_name).variable.ub
        original_lb = tmodel.log_concentration.get_by_id(met_name).variable.lb
        if original_lb < conc_ub:
            tmodel.log_concentration.get_by_id(met_name).variable.ub = conc_ub
        else:
            print(met_name)
            tmodel.log_concentration.get_by_id(met_name).variable.ub = original_lb * 0.99
        try:
            sol = tmodel.optimize()
        except:
            print('cannot impose concentration bound on ' + met_name)
            failed_mets.append(met_name)
            tmodel.log_concentration.get_by_id(met_name).variable.ub = original_ub
            failed_count += 1

# try to minimize the failed mets and also maximize their fluxes
tva_min = []
tva_max = []
for i in failed_mets:
    tva_min.append(_variability_analysis_element(tmodel, tmodel.variables['LC_' + i], 'min'))
    tva_max.append(_variability_analysis_element(tmodel, tmodel.variables['LC_' + i], 'max'))

tva_lc = {}
j = 0
for i in failed_mets:
    tva_lc[i] = [tva_min[j], tva_max[j]]
    j += 1

tmodel.objective = tmodel.reactions.biomass
lc_ub = max(tva_min)
failed_count_2 = 0
failed_mets_2 = []
for met_name in tva_lc.keys():
    original_ub = tmodel.log_concentration.get_by_id(met_name).variable.ub
    tmodel.log_concentration.get_by_id(met_name).variable.ub = lc_ub
    try:
        sol = tmodel.optimize()
    except:
        print('cannot impose concentration bound on ' + met_name)
        failed_mets.append(met_name)
        tmodel.log_concentration.get_by_id(met_name).variable.ub = original_ub
        failed_count += 1