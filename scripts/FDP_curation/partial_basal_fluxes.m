function model = partial_basal_fluxes(model,time_limit,cplex_param,basal_flux, exlude_inds)
%PARTIAL_BASAL_FLUXES imposes basal flux constraint to the specified fluxes
%of the model
%   model: model to impose basal fluxes
%   time_limit: for the solver
%   cplex_param: customize the solver's parameters
%   basal_flux: value for the constraint
%   exlude_inds: ids of fluxes to exlude

[mipTolInt,scalPar,feasTol,emphPar] = setCplexParamHelper(cplex_param);
NF_PNF = getAllVar(model, {'NF', 'PERTURB_NF'});

for i =1:length(model.rxns)
    % Ignore the exlude_inds
    check = find(ismember(exlude_inds, NF_PNF(i)));
    if ~isempty(check)
        continue
    end
    % Impose constraint such that Fi + Ri >= basal_flux
    ind_F=find(ismember(model.varNames,strcat('F_',model.rxns{i})));
    ind_R=find(ismember(model.varNames,strcat('R_',model.rxns{i})));
    [num_constr,~] = size(model.A);
    model.rhs(num_constr+1,1)             = basal_flux;
    model.constraintNames{num_constr+1,1} = strcat('basal_flux',model.rxns{i});
    model.constraintType{num_constr+1,1}  = '>';
    model.A(num_constr+1,[ind_F,ind_R])         = 1;
end

% In the case we have a REMI model
ind_P_F_R = getAllVar(model,{'PERTURB_F';'PERTURB_R'});
if ~isempty(ind_P_F_R)
    for i =1:length(model.rxns)
        % Ignore the exlude_inds
        check = find(ismember(exlude_inds, NF_PNF(length(model.rxns)+i)));
        if ~isempty(check)
            continue
        end
        % Impose constraint such that P_Fi + P_Ri >= basal_flux
        ind_F=find(ismember(model.varNames,strcat('PERTURB_F_',model.rxns{i})));
        ind_R=find(ismember(model.varNames,strcat('PERTURB_R_',model.rxns{i})));
        [num_constr,~] = size(model.A);
        model.rhs(num_constr+1,1)             = basal_flux;
        model.constraintNames{num_constr+1,1} = strcat('perturb_basal_flux',model.rxns{i});
        model.constraintType{num_constr+1,1}  = '>';
        model.A(num_constr+1,[ind_F,ind_R])         = 1;
    end
end

% Check if model is feasible
sol = solveTFAmodelCplex(model,time_limit,[],mipTolInt,emphPar,feasTol,scalPar);
if isempty(sol.x)
    disp('Assigning the basal fluxes resulted in an infeasible model!')
end
end

