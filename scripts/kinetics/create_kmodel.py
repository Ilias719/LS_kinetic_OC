"""
Code to create kinetic model based on the pytfa model provided
The produced kinetic model is a draft version and might be unusable if:
i) The reactant - product pairs are not aligned
ii) the hill_coefficients are set to none and not a value like 1.0

After doing manual changes the model is ready for sampling of kinetic parameters
"""

from skimpy.analysis.oracle.minimum_fluxes import MinFLux, \
    MinFLuxVariable
from pytfa.io.json import load_json_model
from pytfa.analysis import variability_analysis

from skimpy.io.generate_from_pytfa import FromPyTFA
from skimpy.io.yaml import export_to_yaml
from skimpy.core.compartments import Compartment
from skimpy.utils.general import sanitize_cobra_vars

# Thermo model path
path_to_tmodel = '../../models/tfa_model_ATP_modified_with_biomass_final.json'

# Kinetic model output path
path_to_kmodel = '../../models/kinetic_models/kinetic_model_draft.yaml'

# Import and set solver
tmodel = load_json_model(path_to_tmodel)

CPLEX = 'optlang-cplex'
tmodel.solver = CPLEX

tmodel.solver.configuration.tolerances.feasibility = 1e-9
tmodel.solver.configuration.tolerances.optimality = 1e-9
tmodel.solver.configuration.tolerances.integrality = 1e-9
# tmodel.solver.problem.parameters.read.scale = -1
# tmodel.solver.problem.parameters.emphasis.numerical = 1
# tmodel.solver.configuration.presolve = True

# Convert to kinetic model
sol_fdp = tmodel.optimize()

# Metabolites that are not treated as normal products / reactants
small_molecules = ['h_c', 'h_e', 'h_m', 'h_i', 'h_x', 'h_r', 'h_n']

# Metabolites that will not be part of the kinetics of the model
reactants_to_exclude = ['na1_c', 'na1_e']

# Build the kinetic model
model_gen = FromPyTFA(small_molecules=small_molecules,
                      reactants_to_exclude=reactants_to_exclude,
                      max_revesible_deltag_0=50)

kmodel = model_gen.import_model(tmodel,
                                sol_fdp.raw,
                                concentration_scaling_factor=1e6)

# Add and map compartments
for c in tmodel.compartments:
    comp = Compartment(name=c)
    kmodel.add_compartment(comp)

for met in tmodel.metabolites:
    comp = kmodel.compartments[met.compartment]
    kin_met = sanitize_cobra_vars(met.id)
    if kin_met in kmodel.reactants:
        kmodel.reactants[kin_met].compartment = comp
    if kin_met in kmodel.parameters:
        kmodel.parameters[kin_met].compartment = comp

# Add volume parameters
reference_cell_volume = {'cell_volume_c': 1760.,  # micrometersˆ3 from BioNumbers
                         'cell_volume_e': 1760.,
                         'cell_volume_m': 1760.,
                         'cell_volume_x': 1760.,
                         'cell_volume_r': 1760.,
                         'cell_volume_n': 1760.,
                         'cell_volume_i': 1760.,
                         'volume_c': 0.7 * 1760,  # Luby-Phelps K. 2013
                         'volume_e': 1760.,
                         'volume_m': 0.13 * 1760,  # m 13% of cell volume from BioNumbers
                         'volume_x': 0.01 * 1760.,
                         'volume_r': 0.01 * 1760.,
                         'volume_n': 0.01 * 1760.,
                         'volume_i': 0.01 * 1760.
                         }

kmodel.parameters = reference_cell_volume

kmodel.prepare()

# Export the kinetic model
export_to_yaml(kmodel, path_to_kmodel)
