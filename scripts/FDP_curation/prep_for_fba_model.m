function [model_fba] = prep_for_fba_model(model, threshold)
%PREP_FOR_FBA_MODEL prepares the inputs for the production of the fba
%version of the remi model.
%   model: the tfa remi model with NF and PNF bounds from Variability
%   Analysis.
%   threshold: the threshold used when creating the original REMI model

% Initialize the two models: one for NF and one for PERTURB_NF
model1 = model;
model2 = model;

NF = getAllVar(model, {'NF'});
PNF = getAllVar(model, {'PERTURB_NF'});

model1.lb = model.var_lb(NF);
model1.ub = model.var_ub(NF);

model2.lb = model.var_lb(PNF);
model2.ub = model.var_ub(PNF);

% Find the active reaction names and directions
YBplus = getAllVar(model, {'YBplus'});
active_rxns = model.varNames(YBplus(model.var_lb(YBplus)>0.9));

names = strrep(active_rxns, 'YBplus', '');
names = strrep(names, '_F_', '');
names = strrep(names, '_R_', '');

for i = 1:length(active_rxns)
    if contains(active_rxns(i), 'F_')
        dir{i} = 'F';
    else
        dir{i} = 'R';
    end
end

% Find the values of the ratios
cplus = getAllCons(model, {'Cplus'});

j = 1;
for i = 1:length(active_rxns)
    ind = find(ismember(model.constraintNames(cplus), strrep(active_rxns(i), 'YBplus','Cplus')));
    if ~isempty(ind)
        rows = find(model.A(cplus(ind), :));
        ratios(j) = -full(model.A(cplus(ind), rows(1)));
        j = j +1;
    end
end

% Create the model_fba
model_fba = build_FBA_model_remi(model1,model2, names, dir,ratios',threshold);

end

