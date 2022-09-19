% This code will do the following steps:
%       1) Impose ATP/ADP/AMP metabolomics data with 20% relaxation
%       2) Impose basal fluxes to all reactions that is possible
%       3) Impose minimization of the sum of fluxes
%       4) Impose cancer physiology
%       5) Manually curate the directionality of the intacellular reactions
%       6) Make BDRs with assymetric bounds unidirectional
% 
% Each step checks for solver feasibility and does tva for the BDRs left.
% Toumpe I. 2022

clear, clc

% Load the original model with REMI constraints
cplex_loaded = load_cplex;
changeCobraSolver('cplex_direct')
load('../../models/tfa/model_remi.mat')

%-------------------Step 1: ATP energy charge------------------------------

% Metabolomics data for atp/adp/amp (human cell)
T = readtable('../../ranges/ATP_metabolomics.xlsx');
hs = find(ismember(T{:,7}, 'Homo sapiens'));
atp_id = hs(ismember(T{hs,2}, 'ATP'));
adp_id = hs(ismember(T{hs,2}, 'ADP'));
amp_id = hs(ismember(T{hs,2}, 'AMP'));
atp_value = mean(T{atp_id,3});
adp_value = mean(T{adp_id,3});
amp_value = mean(T{amp_id,3});

% Allow for 20% relaxation of the concnentrations as the energy charge is
% between 0.91 and 0.96, (general ec reports: 0.7-0.95).

slack = 0.2;

id_atp = find(ismember(model.varNames, {'LC_atp_c', 'PERTURB_LC_atp_c'}));
id_adp = find(ismember(model.varNames, {'LC_adp_c', 'PERTURB_LC_adp_c'}));
id_amp = find(ismember(model.varNames, {'LC_amp_c', 'PERTURB_LC_amp_c'}));

model.var_lb(id_atp) = log(atp_value*(1-slack));
model.var_ub(id_atp) = log(atp_value*(1+slack));

model.var_lb(id_adp) = log(adp_value*(1-slack));
model.var_ub(id_adp) = log(adp_value*(1+slack));

model.var_lb(id_amp) = log(amp_value*(1-slack));
model.var_ub(id_amp) = log(amp_value*(1+slack));

sol = solveTFAmodelCplex(model);

%-------------------Step 2: Impose basal fluxes----------------------------

% Find all alternatives that cannot carry basal flux
basal_flux = 1e-6;
[model_tmp, prob] = optimize_basal(model, 60, 'LCSBDefault', basal_flux);

% Apply basal fluxes to all reactions except r0122 and PPAm
ids_to_exlude = find(ismember(model.varNames, {'NF_PPAm', 'NF_r0122', 'PERTURB_NF_PPAm', 'PERTURB_NF_r0122'}));
model_basal = partial_basal_fluxes(model, 60, 'LCSBDefault', basal_flux, ids_to_exlude);
sol_basal = solveTFAmodelCplex(model_basal);

% Do variability analysis 
NF_PNF = getAllVar(model_basal, {'NF', 'PERTURB_NF'});
tva_basal = runTMinMax(model_basal, model_basal.varNames(NF_PNF),180);
save('../../ranges/tva_basal.mat', 'tva_basal')

% Find the BDRs of the system and apply the new bounds
BDR_basal = NF_PNF(tva_basal(:,1)<-1e-9 & tva_basal(:,2)>1e-9);

model_basal.var_lb(NF_PNF) = tva_basal(:,1);
model_basal.var_ub(NF_PNF) = tva_basal(:,2);

%---------------Step 3: Impose minimization of sum of fluxes---------------

model_basal.var_lb(find(model_basal.f)) = sol_basal.x(find(model_basal.f))*0.8;
model_min_sum = min_sum_fluxes(model_basal, 100,'LCSBDefault');

% Do variability analysis only for the previous BDRs
tva_min_sum_fluxes = runTMinMax(model_min_sum, model_min_sum.varNames(BDR_basal),180);
save('../../ranges/tva_min_sum_of_fluxes.mat', 'tva_min_sum_fluxes')

% Find remaining BDRs and apply the new bounds
BDR_min_sum = BDR_basal(tva_min_sum_fluxes(:,1)<-1e-9 & tva_min_sum_fluxes(:,2)>1e-9);
uni_min_sum = setxor(BDR_basal, BDR_min_sum);

model_min_sum.var_lb(uni_min_sum) = tva_min_sum_fluxes(ismember(BDR_basal, uni_min_sum), 1);
model_min_sum.var_ub(uni_min_sum) = tva_min_sum_fluxes(ismember(BDR_basal, uni_min_sum), 2);

%--------------------Step 4: Impose cancer physiology----------------------

% Ensure the TCA cycle goes in one direction
model_ph = model_min_sum;
model_ph.var_lb(find(ismember(model_ph.varNames, 'NF_ACONTm'))) = 1e-6;
model_ph.var_lb(find(ismember(model_ph.varNames, 'PERTURB_NF_ACONTm'))) = 1e-6;
model_ph.var_lb(find(ismember(model_ph.varNames, 'NF_AKGDm'))) = 1e-6;
model_ph.var_lb(find(ismember(model_ph.varNames, 'PERTURB_NF_AKGDm'))) = 1e-6;
model_ph.var_ub(find(ismember(model_ph.varNames, 'NF_SUCOASm'))) = -1e-6;
model_ph.var_ub(find(ismember(model_ph.varNames, 'PERTURB_NF_SUCOASm'))) = -1e-6;
model_ph.var_lb(find(ismember(model_ph.varNames, 'NF_SUCD1m'))) = 1e-6;
model_ph.var_lb(find(ismember(model_ph.varNames, 'PERTURB_NF_SUCD1m'))) = 1e-6;
model_ph.var_lb(find(ismember(model_ph.varNames, 'NF_FUMm'))) = 1e-6;
model_ph.var_lb(find(ismember(model_ph.varNames, 'PERTURB_NF_FUMm'))) = 1e-6;
model_ph.var_lb(find(ismember(model_ph.varNames, 'NF_MDHm'))) = 1e-6;
model_ph.var_lb(find(ismember(model_ph.varNames, 'PERTURB_NF_MDHm'))) = 1e-6;

% Impose physiology and do variability analysis 
[results, tva_physiology] = impose_cancer_physiology(model_ph, 1, BDR_min_sum);
save('../../ranges/tva_physiology.mat', 'tva_physiology')

% Find remaining BDRs and apply the new bounds
BDR_physiology = BDR_min_sum(tva_physiology(:,1)<-1e-9 & tva_physiology(:,2)>1e-9);
uni_physiology = setxor(BDR_min_sum, BDR_physiology);

model_ph.var_lb(uni_physiology) = tva_physiology(ismember(BDR_min_sum, uni_physiology), 1);
model_ph.var_ub(uni_physiology) = tva_physiology(ismember(BDR_min_sum, uni_physiology), 2);

%-----Step 5: Curate the directionality of the intacellular reactions------

% Find the non transport reactions
non_tr_id = find(model.isTrans==0);
NF = getAllVar(model, {'NF'});
non_tr_NF = NF(non_tr_id);

% Find the ids of those non transport reactions
BDR.id = intersect(non_tr_NF,BDR_physiology);
BDR.flux = model.varNames(BDR.id); 
BDR.subsystem = model.subSystems(find(ismember(model.rxns,strrep(BDR.flux,'NF_',''))));

% Find the formula of the BDR
BDR.formula = printRxnFormula(model, replace(model.varNames(BDR.id),'NF_',''));

% Rank by subsystem 
[BDR.subsystem, positions] = sortrows(BDR.subsystem);
BDR.id = BDR.id(positions);
BDR.flux = BDR.flux(positions);
BDR.formula = BDR.formula(positions);

% Save to csv
writetable(struct2table(BDR),'BDRS_basal_draft.csv')

% Assign new directionalities manually and iteratively
T = readtable('BDRS_basal_curated.csv');

model_cur = model_ph;

changed = [];
for i = 1:length(BDR.id)
    disp(BDR.flux(i))
    if ismember(T{i,5},'Forward')
        model_cur.var_lb(BDR.id(i)) = 1e-6;
        sol = solveTFAmodelCplex(model_cur);
        if isempty(sol.x) %if infeasible --> change directionality
            model_cur.var_lb(BDR.id(i)) = model_ph.var_lb(BDR.id(i));
            model_cur.var_ub(BDR.id(i)) = -1e-6;
            sol = solveTFAmodelCplex(model_cur);
            changed = [changed, i];
            if isempty(sol.x)
                disp('Error when imposing flux bound')
                keyboard
            end
        end
    elseif ismember(T{i,5}, 'Reverse')
        model_cur.var_ub(BDR.id(i)) = -1e-6;
        sol = solveTFAmodelCplex(model_cur);
        if isempty(sol.x) %if infeasible --> change directionality
            model_cur.var_ub(BDR.id(i)) = model_ph.var_ub(BDR.id(i));
            model_cur.var_lb(BDR.id(i)) = 1e-6;
            sol = solveTFAmodelCplex(model_cur);
            changed = [changed, i];
            if isempty(sol.x)
                disp('Error when imposing flux bound')
                keyboard
            end
        end
    end
end

% Same process for the perturbed model
PNF = getAllVar(model, {'PERTURB_NF'});

% Find the non transport reactions
non_tr_PNF = PNF(non_tr_id);

% Find the ids of those non transport reactions
P_BDR.id = intersect(non_tr_PNF,BDR_physiology);
P_BDR.flux = model.varNames(P_BDR.id); 
P_BDR.subsystem = model.subSystems(find(ismember(model.rxns,strrep(P_BDR.flux,'PERTURB_NF_',''))));

% Find the formula of the BDR
P_BDR.formula = printRxnFormula(model, replace(model.varNames(P_BDR.id),'PERTURB_NF_',''));

% Rank by subsystem 
[P_BDR.subsystem, P_positions] = sortrows(P_BDR.subsystem);
P_BDR.id = P_BDR.id(P_positions);
P_BDR.flux = P_BDR.flux(P_positions);
P_BDR.formula = P_BDR.formula(P_positions);

% Save to csv
writetable(struct2table(P_BDR),'PERTURB_BDRS_draft.csv')

% Assign new directionalities manually and iteratively
T = readtable('PERTURB_BDRS_curated.csv');

for i = 1:length(P_BDR.id)
    disp(P_BDR.flux(i))
    if ismember(T{i,5},'Forward')
        model_cur.var_lb(P_BDR.id(i)) = 1e-6;
        sol = solveTFAmodelCplex(model_cur);
        if isempty(sol.x) %if infeasible --> change directionality
            model_cur.var_lb(P_BDR.id(i)) = model_ph.var_lb(P_BDR.id(i));
            model_cur.var_ub(P_BDR.id(i)) = -1e-6;
            sol = solveTFAmodelCplex(model_cur);
            changed = [changed, i];
            if isempty(sol.x)
                disp('fatal error')
                keyboard
            end
        end
    elseif ismember(T{i,5}, 'Reverse')
        model_cur.var_ub(P_BDR.id(i)) = -1e-6;
        sol = solveTFAmodelCplex(model_cur);
        if isempty(sol.x) %if infeasible --> change directionality
            model_cur.var_ub(P_BDR.id(i)) = model_ph.var_ub(P_BDR.id(i));
            model_cur.var_lb(P_BDR.id(i)) = 1e-6;
            sol = solveTFAmodelCplex(model_cur);
            changed = [changed, i];
            if isempty(sol.x)
                disp('fatal error')
                keyboard
            end
        end
    end
end

% Do variability analysis 
tva_cur = runTMinMax(model_cur, model_cur.varNames(BDR_physiology), 180);
save('../../ranges/tva_curated.mat', 'tva_cur')

% Find remaining BDRs and apply the new bounds
BDR_cur = BDR_physiology(tva_cur(:,1)<-1e-9 & tva_cur(:,2)>1e-9);
uni_cur = setxor(BDR_physiology, BDR_cur);

model_cur.var_lb(uni_cur) = tva_cur(ismember(BDR_physiology, uni_cur), 1);
model_cur.var_ub(uni_cur) = tva_cur(ismember(BDR_physiology, uni_cur), 2);

%---------------Step 6: Make asymmetric BDRs unidirectional--------------
unbalanced_id = [];
tolerance = 1e-4;
for i = 1:length(BDR_cur)
    lbd = abs(model_cur.var_lb(BDR_cur(i)))-1e-6;
    ubd = abs(model_cur.var_ub(BDR_cur(i)))-1e-6;
    if lbd/ubd < tolerance
        fprintf('Upper bound: %.3e\t Lower bound: %.3e\n',model_cur.var_ub(BDR_cur(i)),model_cur.var_lb(BDR_cur(i)))
        model_cur.var_lb(BDR_cur(i)) = 1e-6;
        unbalanced_id = [unbalanced_id, BDR_cur(i)];
    elseif ubd/lbd < tolerance
        fprintf('Upper bound: %.3e\t Lower bound: %.3e\n',model_cur.var_ub(BDR_cur(i)),model_cur.var_lb(BDR_cur(i)))
        model_cur.var_ub(BDR_cur(i)) = -1e-6;
        unbalanced_id = [unbalanced_id, BDR_cur(i)];
    end
end
sol = solveTFAmodelCplex(model_cur);

% Do final variability analysis for all the fluxes and apply new bounds
tva_final = runTMinMax(model_cur, model_cur.varNames(NF_PNF), 180);
save('../../ranges/tva_final.mat', 'tva_final')

model_final = model_cur;
model_final.var_lb(NF_PNF) = tva_final(:,1);
model_final.var_ub(NF_PNF) = tva_final(:,2);

save('../../models/tfa/model_final.mat', 'model_final')

