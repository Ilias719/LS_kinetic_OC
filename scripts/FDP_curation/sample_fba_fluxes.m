% This code will do the following steps:
%       1) Create an LP version of the model (REMI constraints in FBA)
%       2) Create a set of flux samples with the use of a tweaked ACHR
%       sampler that blocks basal bound areas
%       3) Calcluate the average sample point and assign it's FDP to the
%       remaining BDRs in the model
%
% For more information about the sampler check the function used in the
% script.
% Toumpe I. 2022

clear, clc

cplex_loaded = load_cplex;
changeCobraSolver('cplex_direct')

%----------Load model and produce FBA version of the model-----------------
load('./final_results_atp_modified/model_final.mat')

model_remi = prep_for_fba_model(model_final, 1.3);

% Settings of the sampler
nWarmupPoints = 10000;
nFiles = 10;
par_pools = 16;
nPointsPerFile = 6250;
nStepsPerPoint = 500;
NvZeroTolerance = 1e-8;

[nMet,nRxn] = size(model_remi.S);
fprintf('Original model: %d rxns %d metabolites\n',nRxn,nMet);
modelSampling = model_remi;

%-----------------ACHR sampling that ckecks for basal bounds---------------

fprintf('Create warmup points\n');
% Create warmup points for sampler
warmupPts= createHRWarmup_basal_only(modelSampling, model_final, nWarmupPoints,NvZeroTolerance);
save('../../samples/Warmup_Points' ,'warmupPts')

fprintf('Run sampler for a total of %d steps\n',nFiles*nPointsPerFile*nStepsPerPoint*par_pools);
% Sample model
parfor i = 1:par_pools
    changeCobraSolver('cplex_direct'); % Initialize solver for each pool
    ACHRSampler_basal_only(modelSampling,warmupPts,sprintf('../../samples/basal_sample_%d',i),nFiles,nPointsPerFile,nStepsPerPoint,NvZeroTolerance);
end

% Load samples and assign final FDP
disp('Loading ACHR samples')
batch_points = [];
for i = 1:par_pools
    for j = 1:nFiles
        disp(i)
        try
            load(sprintf('../../samples/1M_samples_atp/basal_sample_%d_%d',i,j))
            batch_points = [batch_points, points];
        catch
            keyboard
        end
    end
end

% Mean of samples
ms = mean(batch_points,2);

NF_PNF = [getAllVar(model_final,{'NF'}); getAllVar(model_final,{'PERTURB_NF'})];
id_BD = find(model_final.var_lb(NF_PNF)<-1e-9 & model_final.var_ub(NF_PNF)>1e-9);

% Assign directionalities to the BDRs
pos_dir = id_BD(ms(id_BD)>0);
neg_dir = id_BD(ms(id_BD)<0);
model_final.var_lb(NF_PNF(pos_dir)) = 1e-6;
model_final.var_ub(NF_PNF(neg_dir)) = -1e-6;
sol = solveTFAmodelCplex(model_final);
if ~sempty(sol.x)
    disp('The mean sample point does not produce a feasible FDP. Check the FDP of its closest sample point!')
else
    model_fdp = model_final;
    save('../../models/tfa/model_final_fdp.mat', 'model_fdp')
end

% At this point use runTminmax to assign the final bounds before moving to
% a single model at the kinetic analysis part of the workflow