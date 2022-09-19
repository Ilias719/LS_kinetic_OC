function [results,tva] = impose_cancer_physiology(model, run_tva, id_BDRs)
%IMPOSE_CANCER_PHYSIOLOGY imposes physiologically relevant optimzation
%targtes and checks if any BDRS became unidirectional
%   model: model to impose physiologically relevant optimzation targets
%   run_tva: set it to 1 for variability analysis
%   id_BDRs: ids for the fluxes to do variability analysis 

if nargin == 1
    run_tva = 0;
    id_BDRs = [];
    tva = nan;
end

% Locate all the targeted reactions
reactions_to_look = {'F_biomass';     'PERTURB_F_biomass'
                     'NF_EX_lac_L_e'; 'PERTURB_NF_EX_lac_L_e'; ...
                     'NF_EX_gln_L_e'; 'PERTURB_NF_EX_gln_L_e'; ...
                     'NF_PGK';        'PERTURB_NF_PGK';        ...
                     'NF_PYK';        'PERTURB_NF_PYK'
                     'NF_ATPS4mi';    'PERTURB_NF_ATPS4mi'};
ind = find(ismember(model.varNames, reactions_to_look));

% Part a: only max growth
sol_a = solveTFAmodelCplex(model);

results = cell(length(reactions_to_look),3);
results(:,1) = model.varNames(ind);

for i=1:length(ind)
    results{i,2} = sol_a.x(ind(i));
end

model.var_lb(find(model.f)) = sol_a.x(find(model.f))*0.8;

% Part b: max ATP production, lactate production, glutamine uptake
model.f(find(model.f)) = 0;

model.f(find(ismember(model.varNames, reactions_to_look(7:8)))) = -1; %max ATP production
model.f(find(ismember(model.varNames, reactions_to_look(9:12)))) = 1; %max ATP production
model.f(find(ismember(model.varNames, reactions_to_look(3:4)))) = 1; %max lactate production
model.f(find(ismember(model.varNames, reactions_to_look(5:6)))) = -1; %max glutamine uptake
sol_b = solveTFAmodelCplex(model);

for i=1:length(ind)
    results{i,3} = sol_b.x(ind(i));
end

% Make sure the solution is not too close to the bounds otherwise the model
% will get too constrained
opt = find(model.f);
for i = 1:length(opt)
    obj = model.f(opt(i));
    if obj==1 && model.var_lb(opt(i))<sol_b.x(opt(i))*0.8
        model.var_lb(opt(i)) = sol_b.x(opt(i))*0.8;
    elseif obj==-1 && model.var_ub(opt(i))>sol_b.x(opt(i))*0.8
        model.var_ub(opt(i)) = sol_b.x(opt(i))*0.8;
    end
end

model.f(find(model.f)) = 0;
model.f(find(ismember(model.varNames, reactions_to_look(1:2)))) = 1;

if run_tva == 1
    tva = runTMinMax(model,model.varNames(id_BDRs),120);
end
end

