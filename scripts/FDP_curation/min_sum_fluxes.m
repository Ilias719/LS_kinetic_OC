function model = min_sum_fluxes(model,time_limit,cplex_param)
%MIN_SUM_FLUXES imposes minimization of the sum of the fluxes constraint to
%the model
%   model: model to impose min_sum_of_fluxes
%   time_limit: for the solver
%   cplex_param: customize the solver's parameters

[mipTolInt,scalPar,feasTol,emphPar] = setCplexParamHelper(cplex_param);

% Keep original optimization criteria intact
sol = solveTFAmodelCplex(model,time_limit,[],mipTolInt,emphPar,feasTol,scalPar);
orig_varLB = model.var_lb(find(model.f));
model.var_lb(find(model.f)) = sol.x(find(model.f));

% Make new objective function that minimizes all fluxes
model_min = model;
model_min.f = zeros(length(model_min.varNames),1);
ind_F_R = getAllVar(model_min,{'F';'R'});
model_min.f(ind_F_R) = 1;
model_min.objtype = 1;
min_flux = solveTFAmodelCplex(model_min,time_limit,[],mipTolInt,emphPar,feasTol,scalPar);
store = min_flux.val;

% In the case we have a REMI model
ind_P_F_R = getAllVar(model_min,{'PERTURB_F';'PERTURB_R'});
if ~isempty(ind_P_F_R)
    model_min.f = zeros(length(model_min.varNames),1);
    model_min.f(ind_P_F_R) = 1;
    model_min.objtype = 1;
    min_flux = solveTFAmodelCplex(model_min,time_limit,[],mipTolInt,emphPar,feasTol,scalPar);
    % Assign constraint sum(perturbed_fluxes) <= 1.01 min_flux
    [num_constr,~] = size(model.A);
    model.rhs(num_constr+1,1)             = 1.01*(min_flux.val);
    model.constraintNames{num_constr+1,1} = 'min_sum_fluxes_P';
    model.constraintType{num_constr+1,1}  = '<';
    model.A(num_constr+1,ind_P_F_R)         = 1;
end


[num_constr,~] = size(model.A);
% Assign constraint sum(fluxes) <= 1.01 min_flux
model.rhs(num_constr+1,1)             = 1.01*(store);
model.constraintNames{num_constr+1,1} = 'min_sum_fluxes';
model.constraintType{num_constr+1,1}  = '<';
model.A(num_constr+1,ind_F_R)         = 1;

% Check if model is feasible and default back to original bounds 
sol = solveTFAmodelCplex(model,time_limit,[],mipTolInt,emphPar,feasTol,scalPar);
if isempty(sol.x)
    disp('Assigning min_sum_of_fluxes resulted in infeasible model!')
end
model.var_lb(find(model.f)) = orig_varLB ;
end
