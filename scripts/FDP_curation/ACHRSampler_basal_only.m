function ACHRSampler_basal_only(model,warmupPoints,fileName,nFiles,pointsPerFile,stepsPerPoint,Nv_tol,initPoint,fileBaseNo,maxTime, basal_tolerance)
% ACHRSAMPLER_BASAL_ONLY Artificial Centering Hit-and-Run sampler with
% basal bound constraints
%
%   model: LP version of the model
%   warmupPoints: warmup points from createHRWarmup_basal_only
%   filename: path to save the samples
%   nFiles: split the samples in files for better post-processing
%   pointsPerFile: the number of points saved per file
%   stepsPerPoint: number of points dicarded bedore every saved point
%   Nv_tol: Satisfy mass balance tolerance
%   initPoint: Initial point (Default = centerPoint)
%   fileBaseNo: Base file number for continuing previous sampler run
%   maxTime: Maximum time limit
%   basal_tolerance: the basal flux constraint imposed
%
% Original script from:
% Markus Herrgard, Gregory Hannum, Ines Thiele, Nathan Price 4/14/06
% 
% ----------------------------MODIFICATIONS--------------------------------
%
% The sampler will also understand if the line it walks on is intersecting
% basal bounds and will split the line in segments. It will then pick based
% on a random chance one of those segments in order to do the
% random walk on
%
%Toumpe I. 2022

warning off MATLAB:divideByZero;

% Set the minimum allowed Nv=0 Tolerance
if (nargin < 7)
    Nv_tol = 1e-8;
end
    
if (nargin < 9)
  fileBaseNo = 0;
else
    if (isempty(fileBaseNo))
        fileBaseNo = 0;
    end
end
if (nargin < 10)
    maxTime = 100*36000;
end
if (nargin < 11)
    basal_tolerance = 1e-6;
end

% Minimum allowed distance to the closest constraint
maxMinTol = 1e-9;
% Ignore directions where u is really small
uTol = 1e-9; 
% Project out of directions that are too close to the boundary
dTol = 1e-14;

% Number of warmup points
[nRxns,nWrmup] = size(warmupPoints);

% Find the center of the space
centerPoint = mean(warmupPoints,2);

%------Not useful if the centerPoint is already outside basal bounds.------
% active_rxns = find(model.lb~=0 & model.ub~=0);
% small = find(1e-6 - abs(centerPoint(active_rxns))>1e-9);
% t_min = 100; %random big number, will be replaced by for loop
% If the centerPoint is too close basal bounds find the closest warmupPoint
% and move towards it until the centerpOint is outside the basal bounds.
% if ~isempty(small)
%     basal_bound = basal_tolerance*sign(centerPoint(small));
%     
%     for i = 1:size(warmupPoints,2)
%         t = (basal_bound-centerPoint(small))./(warmupPoints(small,i)-centerPoint(small));
%         if max(abs(t))<abs(t_min)
%             i_stored = i;
%             t_min = max(t);
%         end
%     end
%     centerPoint = centerPoint + t_min.*(warmupPoints(:,i_stored)-centerPoint);
% end

% Set the start point
if (nargin < 8)
    prevPoint = centerPoint;
else
    if (~isempty(initPoint))
        prevPoint = initPoint;
    else
        prevPoint = centerPoint;
    end
end

fidErr = fopen('ACHRerror.txt','w');

totalStepCount = 0;
totalAcceptStepCount = 0;

h = waitbar(0,'ACHR sampling in progress ...');
totalCount = nFiles*pointsPerFile*stepsPerPoint;

t0 = cputime;
fprintf('File #\tPoint #\tStep #\tTime\t#Time left\n');
for i = 1:nFiles

    % Allocate memory for all points
    points = zeros(nRxns,pointsPerFile); 
    
    pointCount = 1;
    while (pointCount <= pointsPerFile)
            
        % Create the random step size vector
        randVector = rand(stepsPerPoint,1);
        
        stepCount = 1;
        while (stepCount <= stepsPerPoint)
            
            % Pick a random warmup point
            randPointID = ceil(nWrmup*rand);
            randPoint = warmupPoints(:,randPointID);
            
            % Find if the randPoint and the centerPoint belong in opposite
            % planes
            opposite_planes = find(centerPoint.*randPoint<0);
            % Begin the line segmentation based on the number of opposite
            % planes
            if ~isempty(opposite_planes)
                total_points = centerPoint;
                point1 = centerPoint;
                point2 = randPoint;
                while true
                    % Find the closest basal boundary and return the points
                    % that belong to the boundary and are part of the line
                    % too
                    [new_point_close, new_point_far] = find_boundary_points(point1, point2, opposite_planes, basal_tolerance);
                    total_points = [total_points, new_point_close, new_point_far];
                    % The new centerPoint is the point past the basal
                    % boundary
                    opposite_planes = find(new_point_far.*randPoint<0);
                    point1 = new_point_far;
                    % Repeat until the new point is in the same FDP as
                    % the randPoint
                    if isempty(opposite_planes)
                        total_points = [total_points, randPoint];
                        
                        % Have a uniform random chance to pick a pair of
                        % points
                        coin = randsample(1:size(total_points,2)/2,1,true);
                        % Save the centerPoint in case in gets tempoarly
                        % changed
                        centerPoint_original = centerPoint;
                        
                        % in order to not change the original code below,
                        % we set the new points in the place of the
                        % randPoint and centerPoint.
                        centerPoint = total_points(:,coin*2-1);
                        randPoint = total_points(:,coin*2);
                        % Based on line 89
                        if stepCount == 1 && pointCount == 1
                            prevPoint = centerPoint;
                        end
                        break
                   end
                end
            end       
            % Sanity check: make sure that the warmupPoint and the
            % centerPoint are not the same
            if randPoint == centerPoint
                continue
            end
            
            % Get a direction from the center point to the warmup point
            u = (randPoint-centerPoint);
            u = u/norm(u);
            
            % Superficial bounds based on the current FDP the sampler is in
            lb = model.lb;
            lb(sign(model.lb)~=sign(prevPoint)) = basal_tolerance;
            ub = model.ub;
            ub(sign(model.ub)~=sign(prevPoint)) = -basal_tolerance;
            
            % Figure out the distances to upper and lower bounds
            distUb = (ub - prevPoint);
            distLb = (prevPoint - lb);
            
            % Figure out if we are too close to a boundary
            validDir = ((distUb > dTol) & (distLb > dTol));
            %model.rxns(~validDir)
            
            % Figure out positive and negative directions
            posDirn = find(u(validDir) > uTol);
            negDirn = find(u(validDir) < -uTol);
            
            % Figure out all the possible maximum and minimum step sizes
            maxStepTemp = distUb(validDir)./u(validDir);
            minStepTemp = -distLb(validDir)./u(validDir);
            maxStepVec = [maxStepTemp(posDirn);minStepTemp(negDirn)];
            minStepVec = [minStepTemp(posDirn);maxStepTemp(negDirn)];
            
            % Figure out the true max & min step sizes
            maxStep = min(maxStepVec);
            minStep = max(minStepVec);
            %fprintf('%f\t%f\n',minStep,maxStep);
            
            % Find new direction if we're getting too close to a constraint
            if (abs(minStep) < maxMinTol & abs(maxStep) < maxMinTol) | (minStep > maxStep)
                fprintf('Warning %f %f\n',minStep,maxStep);
                continue;
            end
            
            % Pick a rand out of list_of_rands and use it to get a random
            % step distance
            stepDist = randVector(stepCount)*(maxStep-minStep)+minStep;
            
            % Advance to the next point
            try
                curPoint_tmp = prevPoint + stepDist*u;
            catch
                continue
            end
            
            % Make sure the point is viable otherwise discard it
            if basal_tolerance-min(abs(curPoint_tmp(active_rxns)))>1e-9
                continue
            end
            
            ind_eq = find(ismember(model.csense,'E'));
            deviation_ACHR = max(abs(model.S(ind_eq,:)*curPoint_tmp-model.b(ind_eq)));
            
            % If the center point has been temporarly moved, reimpose it
            if exist('centerPoint_original', 'var')
                centerPoint = centerPoint_original;
            end
            
            % Print out errors
            if (mod(totalStepCount,2000)==0)
              fprintf(fidErr,'%10.8f\t%10.8f\t',max(curPoint_tmp-model.ub),max(model.lb-curPoint_tmp));
            end

            timeElapsed = cputime-t0;
            
            % Print step information (not customized for parallel
            % computation)
            if (mod(totalStepCount,5000)==0)  
              timePerStep = timeElapsed/totalStepCount;
              fprintf('%d\t%d\t%d\t%8.2f\t%8.2f\n',i,pointCount,totalStepCount,timeElapsed/60,(totalCount-totalStepCount)*timePerStep/60);
            end
            
            overInd = find(curPoint_tmp > model.ub);
            underInd = find(curPoint_tmp < model.lb);
            
            if (any((model.ub-curPoint_tmp) < 0) || any((curPoint_tmp-model.lb) < 0))
              curPoint_tmp(overInd) = model.ub(overInd);
              curPoint_tmp(underInd) = model.lb(underInd);
            end
            
            if (mod(totalStepCount,2000)==0)
              fprintf(fidErr,'%10.8f\n',full(max(max(abs(model.S*curPoint_tmp-model.b)))));
            end
            
            % Check if curPoint_tmp satisfied Nv_tol
            if (mod(stepCount,10)==0) || (mod(stepCount,stepsPerPoint)==0)
                if (deviation_ACHR < Nv_tol)
                    curPoint = curPoint_tmp;
                    prevPoint = curPoint;
                    stepCount = stepCount + 1;
                else                    
                    randPointID = ceil((nWrmup-1)*rand);
                    prevPoint   = warmupPoints(:,randPointID);
                end
            else                
                prevPoint = curPoint_tmp;
                stepCount = stepCount + 1;
            end
            
            
            % Count the total number of steps
            totalStepCount = totalStepCount + 1;
            if mod(totalStepCount,50) == 0
                waitbar(totalStepCount/totalCount,h);
            end
            
            % Recalculate the center point ONLY if the current point meets Nv_tol
            % need to somehow add this back to the code structure
            if (deviation_ACHR < Nv_tol)
                totalAcceptStepCount = totalAcceptStepCount + 1;            
                centerPoint = ((nWrmup+totalAcceptStepCount)*centerPoint + curPoint_tmp)/(nWrmup+totalAcceptStepCount+1);
            end
            
            % Exit if time exceeded
            if (timeElapsed > maxTime)
                points(:,pointCount) = curPoint;
                file = [fileName '_' num2str(fileBaseNo+i) '.mat'];
                save (file,'points');
%                 save ACHR_last_point.mat curPoint
                close(h);
                return;
            end
            
        end % Steps per point
        
        
        % Add the current point to points
        points(:,pointCount) = curPoint;
       
        pointCount = pointCount + 1;
         
    end % Points per cycle
    
    % Save current points to a file
    file = [fileName '_' num2str(fileBaseNo+i) '.mat'];
    save (file,'points');
    
end
close(h);

% Save last point for restart
% save ACHR_last_point.mat curPoint

fclose(fidErr);