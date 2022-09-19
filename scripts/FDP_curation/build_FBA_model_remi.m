function model_remi = build_FBA_model_remi(model1,model2,rxns_remi,rxns_remi_dir,ratio_rxns_remi,threshold)
%BUILD_FBA_MODEL_REMI creates an FBA version of the REMI model with flux
%ratio constraints
%   model1: model with NF flux bounds
%   model2: model with PNF flux bounds
%   rxns_remi: reaction names with active REMI constraints
%   rxns_remi: reaction directionality with active REMI constraints
%   ratio_rxns_remi: flux ratio based on the REMI constraint
%   threshold: the threshold used when creating the original REMI model

% Build FBA model
model_remi.mets = [strcat('cond1_',model1.mets);strcat('cond2_',model2.mets)];
model_remi.rxns = [strcat('cond1_',model1.rxns);strcat('cond2_',model2.rxns)];
model_remi.S = sparse(length(model_remi.mets),length(model_remi.rxns));
model_remi.S(1:length(model1.mets),1:length(model1.rxns)) = model1.S;
model_remi.S(length(model1.mets)+1:length(model_remi.mets),length(model1.rxns)+1:length(model_remi.rxns)) = model2.S;
model_remi.lb = [model1.lb;model2.lb];
model_remi.ub = [model1.ub;model2.ub];
model_remi.b = [model1.b;model2.b];
model_remi.c = [model1.c;model2.c];
model_remi.csense = repmat('E',length(model_remi.mets),1);

% Add constraint ratios consistent with REMI solution
for i = 1:length(rxns_remi)

    if ratio_rxns_remi(i)>=threshold
        % v_cond2/v_cond1>=ratio; ratio*v_cond1 - v_cond2 <= 0
        [~,ind_v_cond1] = ismember(strcat('cond1_',rxns_remi(i)), model_remi.rxns);
        [~,ind_v_cond2] = ismember(strcat('cond2_',rxns_remi(i)), model_remi.rxns);
        model_remi.mets(end+1) = {strcat('const_ratio_',rxns_remi{i})};
        model_remi.S(end+1,ind_v_cond2) = 1;
        model_remi.S(end,ind_v_cond1) = -ratio_rxns_remi(i);
        model_remi.b(end+1) = -1e-7;
        model_remi.csense(end+1) = 'G';

        if strcmp(rxns_remi_dir{i},'F')
            if model_remi.lb(ind_v_cond1) <= 1e-4
                model_remi.lb(ind_v_cond1) = 1e-4;
            end
            if model_remi.lb(ind_v_cond2) <= 1e-4
                model_remi.lb(ind_v_cond2) = 1e-4;
            end
        elseif strcmp(rxns_remi_dir{i},'R')
            if model_remi.ub(ind_v_cond1) >= -1e-4
                model_remi.ub(ind_v_cond1) = -1e-4;
            end
            if model_remi.ub(ind_v_cond2) >= -1e-4
                model_remi.ub(ind_v_cond2) = -1e-4;
            end
        end
    elseif ratio_rxns_remi(i)<=1/threshold
        % v_cond2/v_cond1<=ratio; v_cond2 - ratio*v_cond1 <= 0
        [~,ind_v_cond1] = ismember(strcat('cond1_',rxns_remi(i)), model_remi.rxns);
        [~,ind_v_cond2] = ismember(strcat('cond2_',rxns_remi(i)), model_remi.rxns);
        model_remi.mets(end+1) = {strcat('const_ratio_',rxns_remi{i})};
        model_remi.S(end+1,ind_v_cond2) = 1;
        model_remi.S(end,ind_v_cond1) = -ratio_rxns_remi(i);
        model_remi.b(end+1) = 1e-7;
        model_remi.csense(end+1) = 'L';
        
        if strcmp(rxns_remi_dir{i},'F')
            if model_remi.lb(ind_v_cond1) <= 1e-4
                model_remi.lb(ind_v_cond1) = 1e-4;
            end
            if model_remi.lb(ind_v_cond2) <= 1e-4
                model_remi.lb(ind_v_cond2) = 1e-4;
            end
        elseif strcmp(rxns_remi_dir{i},'R')
            if model_remi.ub(ind_v_cond1) >= -1e-4
                model_remi.ub(ind_v_cond1) = -1e-4;
            end
            if model_remi.ub(ind_v_cond2) >= -1e-4
                model_remi.ub(ind_v_cond2) = -1e-4;
            end
        end
    else
        disp('something is wrong with the ratios') 
    end

    sol.f = 0;
    try
        sol = solveFBAmodelCplex(model_remi,-1,1e-9);
         while sol.f==0 || isempty(sol.f)
             i
            model_remi.b(end) = 10*model_remi.b(end);
            sol = solveFBAmodelCplex(model_remi,-1,1e-9);
        end
    catch
        while sol.f==0 || isempty(sol.f)
            model_remi.b(end) = 10*model_remi.b(end);
            sol.f = 0;
            try
                sol = solveFBAmodelCplex(model_remi,-1,1e-9);
            catch
                [i model_remi.b(end)]
            end
        end
    end

end
sol = solveFBAmodelCplex(model_remi,-1,1e-9);
 

[~,ind_v_cond1] = ismember(strcat('cond1_',rxns_remi), model_remi.rxns);
[~,ind_v_cond2] = ismember(strcat('cond2_',rxns_remi), model_remi.rxns);
[ratio_rxns_remi  sol.x(ind_v_cond2)./sol.x(ind_v_cond1)]


end


  