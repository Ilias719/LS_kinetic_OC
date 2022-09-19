"""
Code that checks if there are any turnover ratios (sum(consuming fluxes)/(metabolite concentration)) are slower than the
eigenvalue cutoff required for physiologically relevant dynamics. This is done only for the metaboltes in linear
pathways because it can be proved that their eigenvalues are heavily depedent by the respective turnover ratios.
If there are such instances, a relaxation to the connected DGo bounds is performed.
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
from element_variability_analysis import element_variability_analysis

# Load pytfa model
path_to_tmodel = '../../models/tfa/tfa_model_no_interventions.json'
tmodel = load_json_model(path_to_tmodel)
tmodel.solver.configuration.tolerances.feasibility = 1e-9
tmodel.solver.configuration.tolerances.optimality = 1e-9
tmodel.solver.configuration.tolerances.integrality = 1e-9

# Scaling parameters
CONCENTRATION_SCALING = 1e6  # 1 mol to 1 mmol
TIME_SCALING = 1  # 1hour

# Parameters of human cancer cells
DENSITY = 1200  # g/L
GDW_GWW_RATIO = 0.25  # Assumes 75% Water

flux_scaling_factor = 1e-3 * (GDW_GWW_RATIO * DENSITY) * CONCENTRATION_SCALING / TIME_SCALING
eigenvalue_cutoff = 3 / 26.1  # 3 times faster than the doubling time of the cell (26.1 hours)

# Exlude some metabolites that are fixed
metabolites_to_exlude = ['h_c', 'h_e', 'h_m', 'h_i', 'h_x', 'h_r', 'h_n', 'h2o_x']

# Find the metabolites that belong in linear pathways
lin_mets = [i for i in tmodel.metabolites if len(i.reactions) == 2 and i.id not in metabolites_to_exlude]

# Best case scenario:   reaction.upper_bound / metabolite.lower_bound
tc_ratio = {}
bad_mets = []
dg_stds = {}
for ix in lin_mets:
    reac = list(ix.reactions)[0]
    r_ub = max(np.abs(reac.bounds)) * flux_scaling_factor
    conc_lb = np.exp(tmodel.log_concentration.get_by_id(ix.id).variable.lb) * CONCENTRATION_SCALING
    tc_ratio[ix.id] = r_ub / conc_lb
    if r_ub / conc_lb < eigenvalue_cutoff:
        print(ix.id)
        bad_mets.append(ix.id)
        for r in ix.reactions:
            dg_stds[r.id] = element_variability_analysis(tmodel, 'DGo_{}'.format(r.id))

# Relax Dgo bounds for the lin_mets with slow TCs
j = 0
dgo_tol = 5
for i in dg_stds:
    ranges = dg_stds[i]
    if ranges[0] > dgo_tol:
        tmodel.log_concentration.get_by_id(bad_mets[int(np.floor(j/2))]).variable.lb = -26
        tmodel.delta_gstd.get_by_id(i).variable.lb = -dgo_tol
        tmodel.delta_gstd.get_by_id(i).variable.ub = dgo_tol
        tmodel.delta_g.get_by_id(i).variable.lb = -1000
        tmodel.delta_g.get_by_id(i).variable.ub = 1000
    j += 1
