function proj_points = pflat(points)
%------------------------------------------------------------------
% PURPOSE
%
%   Project points from given matrix
%
% INPUT:   points:              Matrix [m x n]
%
% OUTPUT   proj_points:         Matrix [m x n-1]
%------------------------------------------------------------------

%   Nils Broman, 2021-01-22
%------------------------------------------------------------------

%proj_points=points(1:end-1,:)./points(end,:);
proj_points=points./points(end,:);
end