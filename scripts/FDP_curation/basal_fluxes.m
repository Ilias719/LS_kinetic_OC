function model = basal_fluxes(model,time_limit,cplex_param,basal_flux)

[mipTolInt,scalPar,feasTol,emphPar] = setCplexParamHelper(cplex_param);

for i =1:length(model.rxns)
    ind_F=find(ismember(model.varNames,strcat('F_',model.rxns{i})));
    ind_R=find(ismember(model.varNames,strcat('R_',model.rxns{i})));
    [num_constr,~] = size(model.A);
    model.rhs(num_constr+1,1)             = basal_flux;
    model.constraintNames{num_constr+1,1} = strcat('basal_flux',model.rxns{i});
    model.constraintType{num_constr+1,1}  = '>';
    model.A(num_constr+1,[ind_F,ind_R])         = 1;
end

% in the case we have a remi model
ind_P_F_R = getAllVar(model,{'PERTURB_F';'PERTURB_R'});
if ~isempty(ind_P_F_R)
    for i =1:length(model.rxns)
        ind_F=find(ismember(model.varNames,strcat('PERTURB_F_',model.rxns{i})));
        ind_R=find(ismember(model.varNames,strcat('PERTURB_R_',model.rxns{i})));
        [num_constr,~] = size(model.A);
        model.rhs(num_constr+1,1)             = basal_flux;
        model.constraintNames{num_constr+1,1} = strcat('perturb_basal_flux',model.rxns{i});
        model.constraintType{num_constr+1,1}  = '>';
        model.A(num_constr+1,[ind_F,ind_R])         = 1;
    end
end

sol = solveTFAmodelCplex(model,time_limit,[],mipTolInt,emphPar,feasTol,scalPar);
if isempty(sol.x)
    disp('Assigning the basal fluxes resulted in an infeasible model!')
end

end

