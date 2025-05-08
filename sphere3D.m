% This function computes the boundary points and normal vectors for a
% sphere using a Fibonacci lattice
%
% radius = radius of sphere
% points = number of points
%
% Written by J. Simpson on 5/8/2025

function [bdy, normal] = sphere3D(radius,points)

phi = zeros(1,points);
theta = zeros(1,points);
golden_ratio = (1+sqrt(5))/2;

for i = 1:points
    theta(i) = 2*pi*mod((i-1)/golden_ratio,1);
    phi(i) = acos(1-2*(i-0.5)/(points+1));
end

% boundary points

bdy = [radius.*sin(phi).*cos(theta); radius.*sin(phi).*sin(theta); radius.*cos(phi)];

% normal vectors

normal = [sin(phi).*cos(theta); sin(phi).*sin(theta); cos(phi)];

end

