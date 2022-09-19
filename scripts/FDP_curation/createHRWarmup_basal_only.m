function warmupPts = createHRWarmup_basal_only(model, model_tfa, basal_tolerance, nPoints,Nv_tol,verbFlag,nPointsCheck)
% CREATEHRWarmup_basal_only Create a warmup point set for hit-and-run
% sampling by combining orthogonal and random points
%
% Additional modification in order to ensue basal fluxes for the BDR
% reactions:
% Step 1: Find solution with MILP model
% Step 2: Extract the FDP and apply it to the FBA model
% Step 3: Find solution with LP model
% Step 4: Reset the FBA model's bounds
% 
%   model: LP version of the model
%   model_tfa: MILP version of the model
%   basal_tolerance: the basal flux constraint imposed
%   nPoints: Number of warmup points
%   Nv_tol: Satisfy mass balance tolerance
%   verbFlag:  Verbose flag
%
% Original script from Markus Herrgard, Richard Que
%
% Toumpe I. 2022: Added basal fluxes constraint for sampled warmup points

if (nargin < 4)||isempty(nPoints), nPoints = 5000; end
if (nargin < 5)||isempty(Nv_tol), Nv_tol = 1e-8; end
if (nargin < 6)||isempty(verbFlag), verbFlag = true; end
if (nargin < 7)||isempty(nPointsCheck), nPointsCheck = true; end

if isfield(model,'S')
    [nMets,nRxns] = size(model.S);
    model.A=model.S;
elseif isfield(model,'A')
    [nMets,nRxns] = size(model.A);
end
if ~isfield(model,'csense')
    model.csense(1:size(model.S,1)) = 'E';
end

if nPointsCheck && (nPoints < nRxns*2) 
    warning(['Need a minimum of ' num2str(nRxns*2) ' warmup points']);
    nPoints = nRxns*2;
end

% cplex_loaded = load_cplex;
% changeCobraSolver('cplex_direct')

% Initialize variables and parameters
warmupPts = sparse(nRxns,nPoints);
model_orig = model;
i = 1;
i_store = 0;
h = waitbar(0,'Creating warmup points ...');
NF_PNF = getAllVar(model_tfa, {'NF','PERTURB_NF'});
id_BD = find(model_tfa.var_lb(NF_PNF)<-1e-9 & model_tfa.var_ub(NF_PNF)>1e-9);
model_tfa.f(model_tfa.f==1) = 0;

% Find the fluxes that are blocked and ignore them
blocked_rxns = find(model_tfa.var_lb(NF_PNF)==0 & model_tfa.var_ub(NF_PNF)==0);
active_rxns = 1:nRxns;
active_rxns(blocked_rxns) = [];

% Generate the points
while i_store < nPoints
    if mod(i,10) == 0
        waitbar(i_store/nPoints,h);
    end
    
    % Create random objective function
    model.c = rand(nRxns,1)-0.5;
    
    for maxMin = [1, -1]
        % Set the objective function
        if i <= nRxns
            model.c = zeros(nRxns,1);
            model.c(i) = 1;
        end
        model.osense = maxMin;
        
        % Step 1: Solve with the MILP model
        model_tfa.f(NF_PNF) = model.c;
        model_tfa.objtype = maxMin;
        sol_tfa = solveTFAmodelCplex(model_tfa, 30);
        
        if ~isempty(sol_tfa.x)
            % Step 2: Extract the FDP and apply it to the FBA model
            direction = sign(sol_tfa.x(NF_PNF(id_BD)));
            model.lb(id_BD(direction == 1)) = basal_tolerance;
            model.ub(id_BD(direction == -1)) = -basal_tolerance;
            clear sol_tfa
        end
        % Step 3: Solve FBA model
        
        % Determine the max or min for the rxn
        sol = solveCobraLP(model);
        x = sol.full;
        status = sol.stat;
        if status == 1
            validFlag = true;
        else
            disp('invalid solution')
            validFlag = false;
            display(status)
            pause;
        end       
        
        % There might be non-blocked reactions that still do not satisfy
        % basal flux constraint.
        prob = find(basal_tolerance-abs(x(active_rxns))>1e-9);

        if ~isempty(prob) 
            prob = active_rxns(prob);
            while true
                % Extract the directionality of the solution and enforce
                % basal bounds
                model.lb(prob(sign(x(prob))==1)) = basal_tolerance;
                model.ub(prob(sign(x(prob))==-1)) = -basal_tolerance;
                % In the case where the reaction is 0, assign a random
                % directionality
                blocked = prob(sign(x(prob))==0);
                if ~isempty(blocked)
                    for ind = 1:length(blocked)
                        coin = round(rand);
                        if coin == 0
                            model.lb(blocked(ind))= basal_tolerance;
                        else
                            model.ub(blocked(ind))= -basal_tolerance;
                        end
                    end
                end
                % Solve based on the new FDP
                sol = solveCobraLP(model);
                x = sol.full;
                status = sol.stat;
                if status == 1
                    validFlag = true;
                else
                    disp('invalid solution')
                    validFlag = false;
                    display(status)
                    pause;
                end  
                % Check if new fluxes do not satisfy the basal tolerance
                prob = find(basal_tolerance-abs(x)>1e-9, 1);
                if isempty(prob)
                    break
                end
            end
        end
            
        % Continue if optimal solution is found     
        
        % Move points to within bounds
        x(x > model.ub) = model.ub(x > model.ub);
        x(x < model.lb) = model.lb(x < model.lb);
        
        ind_eq = find(ismember(model.csense,'E'));
        deviation_wmpPts = max(abs(model.S(ind_eq,:)*x-model.b(ind_eq)));
        
        % Store point only if N*v Tolerance is satisfied
        if (maxMin == 1) 
            if (deviation_wmpPts < Nv_tol) && validFlag
                i_store = i_store+1;
                warmupPts(:,i_store) = x;
            else
                display ('N*v Tolerance not satisfied')
                fprintf('%f\n',i_store/i/2);
            end    
        else
            if (deviation_wmpPts < Nv_tol) && validFlag && i_store < nPoints   
                i_store = i_store+1;
                warmupPts(:,i_store) = x;            
            else
                display ('N*v Tolerance not satisfied')
                fprintf('%f\n',i_store/i/2);
            end
        end
        
        if (verbFlag)
            if mod(i,100)==0
                fprintf('%4.1f\n',i_store/nPoints*100);
            end
        end
        
        % Step 4: Reset the FBA model's bounds
        model = model_orig; 
    end
    if validFlag
        i = i+1;
    end 
    % This if statement is important to continue with the basal flux
    % constraints on samples but could be implemented better. After the
    % orthogonal sampling, it checks if the centerPoint is inside basal
    % bounds. If it is, it continues to sample. The ACHR_basal_only code
    % requires that the centerPoint is away from the basal bounds.
    if i_store >=2*nRxns
        centerPoint = mean(warmupPts(:,1:i_store),2);
        if basal_tolerance - min(abs(centerPoint(active_rxns))) < 1e-9
            break
        else
            disp('Sampling finished but the centerPoint is inside the basal bounds. Sample more warmup points')
        end
    end
end
warmupPts = warmupPts(:,1:i_store);
centerPoint = mean(warmupPts,2);
% Move points in for the ACHR sampling process
save('../../samples/Warmup_Points_no_center_movement' ,'warmupPts')
disp('Moving warmup points closer to the centerPoint')
warmupPts = full(warmupPts);
active_rxns = find(centerPoint ~=0);
for i = 1:size(warmupPts,2)
    curPoint = warmupPts(:,i);
    guess = curPoint*.33 + .67*centerPoint;
    % If the moved warmup point does not satisfy the basal tolerance, we
    % find the closest boundary point based on the line between the warmup
    % point and the centerpoint.
    if basal_bound-min(abs(guess(active_rxns)))>1e-10
        opposite_points = find(centerPoint.*curPoint<0);
        % Begin the line segmentation based on the number of opposite
        % planes
        if ~isempty(opposite_points)
            total_points = centerPoint;
            point1 = centerPoint;
            point2 = curPoint;
            while true
                % Find the closest basal boundary and return the points
                % that belong to the boundary and are part of the line.
                [new_point_close, new_point_far] = find_boundary_points(point1, point2, opposite_points, basal_tolerance);
                total_points = [total_points, new_point_close, new_point_far];
                % The new centerPoint is the point past the basal
                % boundary.
                opposite_points = find(new_point_far.*curPoint<0);
                point1 = new_point_far;
                % Repeat until the new point is in the same FDP as
                % the randPoint.
                if isempty(opposite_points)
                    total_points = [total_points, curPoint];
                    guess = full(guess);
                    diff = [];
                    for j = 1:size(total_points,2)
                        diff(j) = norm(total_points(:,j)-guess);
                    end
                    % Choose the point which is closest to the guess
                    % variable
                    ind = find(diff == min(diff));
                    ind = ind(1);
                    if ind == 1 || ind == length(diff)
                        disp('the closest point is not acceptable!')
                        keyboard
                    end
                    warmupPts(:,i) = total_points(:,ind);
                    break
                end
            end
        end
    else
        warmupPts(:,i) = guess;
    end
end


close(h);
