% This function computes the boundary points and normal vectors for two
% spheres using a Fibonacci lattice
%
% radius1 = radius of sphere1
% radius2 = radius of sphere2
% center1 = center of sphere1
% center2 = center of sphere2
% points1 = number of points of sphere1
% points2 = number of points of sphere2
%
% Written by J. Simpson on 5/8/2025

function [bdy, normal] = twospheres(radius1,radius2,center1,center2,points1,points2)

phi1 = zeros(1,points1);
phi2 = zeros(1,points2);
theta1 = zeros(1,points1);
theta2 = zeros(1,points2);
golden_ratio = (1+sqrt(5))/2;

for i = 1:points1
    theta1(i) = 2*pi*mod((i-1)/golden_ratio,1);
    phi1(i) = acos(1-2*(i-0.5)/(points1+1));
end

for i = 1:points2
    theta2(i) = 2*pi*mod((i-1)/golden_ratio,1);
    phi2(i) = acos(1-2*(i-0.5)/(points2+1));
end

% boundary points

sphere1 = [radius1.*sin(phi1).*cos(theta1)+center1(1);...
    radius1.*sin(phi1).*sin(theta1)+center1(2);...
    radius1.*cos(phi1)+center1(3)];
sphere2 = [radius2.*sin(phi2).*cos(theta2)+center2(1);...
    radius2.*sin(phi2).*sin(theta2)+center2(2);...
    radius2.*cos(phi2)+center2(3)];

% normals vectors

normal1 = [sin(phi1).*cos(theta1); sin(phi1).*sin(theta1); cos(phi1)];
normal2 = [sin(phi2).*cos(theta2); sin(phi2).*sin(theta2); cos(phi2)];

bdy = [sphere1,sphere2];
normal = [normal1,normal2];