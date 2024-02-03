% This function finds the projection of a vector onto a line passing
% through the origin
% Inputs:
% 1) p1: Base of the vector p1 = [x1, y1]  x1: x-coordinate, y1: y-coordinate
% 2) p2: Head of the vector p2 = [x2, y2]
% 3) slope: Slope of the line which passes through the origin onto which we
%           want to project our vector

% Outputs:
% 1) proj: Projection onto the line. If this value is negative, it is
% projected in the direction of -inf 

function [proj] = FindProjection(vec, slope)
line = [1, slope];
proj = zeros(1,size(vec,2));
for i = 1:size(vec, 2)
    proj(i) = dot(line,vec(:,i))/dot(line, line);
end
end