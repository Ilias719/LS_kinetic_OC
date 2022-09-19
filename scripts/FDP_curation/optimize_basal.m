function [model, problematic_reactions] = optimize_basal(model,time_limit,cplex_param,basal_flux)
%OPTIMIZE_BASAL uses MILP formulation to find the maximum possible imposed
%basal flux constraints as well as alternatives
%   model: model to impose basal fluxes
%   time_limit: for the solver
%   cplex_param: customize the solver's parameters
%   basal_flux: value for the constraint

[mipTolInt,scalPar,feasTol,emphPar] = setCplexParamHelper(cplex_param);

IndOfDebugVars = [];
orig_f = find(model.f);
model.f = zeros(length(model.f),1);
[~,num_vars] = size(model.A);

% Assign binary variables
for i = 1:length(model.rxns)
    model.varNames{num_vars+i,1} = strcat('Binary_',model.rxns{i});
    model.var_ub(num_vars+i,1) = 1;
    model.var_lb(num_vars+i,1) = 0;
    model.vartypes{num_vars+i,1} = 'B';
    % Set the vector of the objective function for the newly generated
    % variables to 1
    model.f(num_vars+i,1)=1; 
    IndOfDebugVars(i,1)=num_vars+i;
end

% Assign the basal flux constraint F + R + basal_flux*B > basal_flux
for i =1:length(model.rxns)
    ind_F=find(ismember(model.varNames,strcat('F_',model.rxns{i})));
    ind_R=find(ismember(model.varNames,strcat('R_',model.rxns{i})));
    [num_constr,~] = size(model.A);
    model.rhs(num_constr+1,1)             = basal_flux;
    model.constraintNames{num_constr+1,1} = strcat('basal_flux_',model.rxns{i});
    model.constraintType{num_constr+1,1}  = '>';
    model.A(num_constr+1,[ind_F,ind_R])         = 1;
    model.A(num_constr+1,IndOfDebugVars(i))         = +basal_flux;
end

% In the case we have a REMI model
ind_P_F_R = getAllVar(model,{'PERTURB_F';'PERTURB_R'});
if ~isempty(ind_P_F_R)
    % Assign binary variables for Perturbed model
    [~,num_vars] = size(model.A);
    for i = 1:length(model.rxns)
        model.varNames{num_vars+i,1} = strcat('PERTURB_Binary_',model.rxns{i});
        model.var_ub(num_vars+i,1) = 1;
        model.var_lb(num_vars+i,1) = 0;
        model.vartypes{num_vars+i,1} = 'B';
        % Set the vector of the objective function for the newly generated
        % variables to 1
        model.f(num_vars+i,1)=1; 
        IndOfDebugVars(end+1,1)=num_vars+i;
    end
    % Assign the basal flux constraint P_F + P_R + basal_flux*P_B > basal_flux
    for i =1:length(model.rxns)
        ind_F=find(ismember(model.varNames,strcat('PERTURB_F_',model.rxns{i})));
        ind_R=find(ismember(model.varNames,strcat('PERTURB_R_',model.rxns{i})));
        [num_constr,~] = size(model.A);
        model.rhs(num_constr+1,1)             = basal_flux;
        model.constraintNames{num_constr+1,1} = strcat('perturb_basal_flux_',model.rxns{i});
        model.constraintType{num_constr+1,1}  = '>';
        model.A(num_constr+1,[ind_F,ind_R])         = 1;
        model.A(num_constr+1,IndOfDebugVars(length(model.rxns)+i))         = +basal_flux;
    end
end

% Solver should minimize the number of the new activated binary variables
model.objtype = 1;

% Solve the system to find the maximum basal fluxes that can be imposed
sol = solveTFAmodelCplex(model,time_limit,[],mipTolInt,emphPar,feasTol,scalPar);

% Get the newly added binary variables that were not able to be set to 0
ValuesofBinaries = sol.x(IndOfDebugVars);

ids_ConstraintsRelaxed = find(ValuesofBinaries==1);
ids_ConstraintsImposed = find(ValuesofBinaries==0);

fprintf('basal fluxes not enforced \n')
disp(model.varNames(IndOfDebugVars(ids_ConstraintsRelaxed)))
problematic_reactions = strrep(model.varNames(IndOfDebugVars(ids_ConstraintsRelaxed)), 'Binary_','');

% Add a constraint that forces the sum(Bi) = sum(ValuesofBinaries)
[num_constr,~] = size(model.A);
model.rhs(num_constr+1,1)             = sum(ValuesofBinaries);
model.constraintNames{num_constr+1,1} = 'maximal binaries enforced';
model.constraintType{num_constr+1,1}  = '=';
model.A(num_constr+1,IndOfDebugVars)         = ones(length(IndOfDebugVars),1);

% Have a save point for the original solution
model_orig = model;
model_orig.objtype = -1;
model_orig.f = zeros(length(model_orig.f),1);
model_orig.f(orig_f) = 1;


% Enumerate alternative solutions
i = 0;
while true
    i = i+1;
    % Add a constraint such that sum(Bi*Bi_alternative) <= sum(Bi_alternative) -1
    [num_constr,~] = size(model.A);
    model.rhs(num_constr+1,1)             = sum(ValuesofBinaries)-1;
    model.constraintNames{num_constr+1,1} = sprintf('alternative solution %d',i);
    model.constraintType{num_constr+1,1}  = '<';
    model.A(num_constr+1,IndOfDebugVars)         = ValuesofBinaries;
    sol = solveTFAmodelCplex(model);
    if isempty(sol.x)
        disp('no more alternatives')
        % Keep only the first solution
        model = model_orig;
        break
    end
    ValuesofBinaries = sol.x(IndOfDebugVars);
    ids_ConstraintsRelaxed = find(ValuesofBinaries==1);
    fprintf('basal fluxes not enforced \n')
    disp(model.varNames(IndOfDebugVars(ids_ConstraintsRelaxed)))
end
end

