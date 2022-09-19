function [new_point_close, new_point_far] = find_boundary_points(centerPoint,randPoint, opposite_planes, basal_tolerance)
%FIND_BOUNDARY_POINTS given two points, this script calculates the line
%between them and finds the closest boundary points
% LINE EQUATION USED: NEW_POINT = A + scalar_value * (B - A)
%
%   centerPoint: starting point
%   randPoint: ending point
%   opposite_planes: where centerPoint.*randPoint < 0
%   basal_tolerance: value of the basal constraint

% Ignore blocked reactions
ar = find(centerPoint~=0); % active points
flag = 0;
basal_bounds = sign(centerPoint).*ones(length(centerPoint),1).*basal_tolerance;
u = (randPoint-centerPoint);

% x_point = centerPoint + t*u
t = (basal_bounds(opposite_planes)-centerPoint(opposite_planes))./u(opposite_planes);
% Find the negative or zero values and snap to the other side of the basal
% bound
t_pos = t(t>0);
other_ind = find(t<=0);
if ~isempty(other_ind)
    basal_bounds(opposite_planes(other_ind)) = basal_tolerance*(-sign(basal_bounds(opposite_planes(other_ind))));
    t_new = (basal_bounds(opposite_planes(other_ind))-centerPoint(opposite_planes(other_ind)))./u(opposite_planes(other_ind));
    t_new = t_new(1);
    new_point_close = centerPoint;
    new_point_far = centerPoint + t_new.*u;
    flag = 1;
end


if flag == 0 %If not in special case
    [vx, ~] = sort(t_pos);
    min_val = vx(1);
    min_pos = find(t == min_val);
    k = 1;
    % Make sure that the point is not one of the blocked reactions
    while true
        if ismember(min_pos,ar)
            break
        else
            k = k + 1;
            min_val = vx(k);
            min_pos = find(t == min_val);
        end
    end
    new_point_close = centerPoint + min_val.*u;

    ind = find(basal_tolerance-abs(new_point_close(ar))>1e-9);
    if ~isempty(ind)
        ind = ar(ind); % global index
        new_point_close(ind) = basal_tolerance.*sign(new_point_close(ind));
    end
    
    % Calculate the basal point from the other side of the closest basal
    % bound
    t_other = (-basal_bounds(opposite_planes(min_pos))-centerPoint(opposite_planes(min_pos)))./u(opposite_planes(min_pos));
    t_other = t_other(1);
    if t_other == 0 % That means that we have flipped a basal bound and we are back at the centerPoint
        new_point_far = new_point_close;
        new_point_close = centerPoint;
    else
        new_point_far = centerPoint + t_other.*u;
        in_bound = [];
        in_bound = find((basal_tolerance-abs(new_point_far(ar)))>1e-9);
        % If the new point lands inside a basal bound find the next closest
        % point to it that is also part of the line
        if ~isempty(in_bound)
           in_bound = ar(in_bound);
           while true
                new_basal_bounds = sign(randPoint).*ones(length(randPoint),1).*basal_tolerance;
                t_new = (new_basal_bounds(in_bound)-centerPoint(in_bound))./(u(in_bound));
                t_new_pos = t_new(t_new>0);
                [vx, ~] = sort(t_new_pos);
                try
                    max_val = vx(end);
                catch
                     max_val = 0;
                end
                new_point_far = centerPoint +max_val.*u;
                in_bound = find((basal_tolerance-abs(new_point_far(ar)))>1e-9);
                in_bound = ar(in_bound);
                if isempty(in_bound)
                    break
                end
           end
        end
    end
end

end

