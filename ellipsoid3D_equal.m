% This function computes the boundary points and normal vectors for the
% ellipsoid satisfying
%       x^2/(1.2^2)+y^2/(1.2)^2+z^2/(0.7^2)=1
% using roughly equally spaced points
%
% scale = relative size of ellipsoid
% level = number of discretization points in z
%
% Written by J. Simpson on 5/8/2025

function [bdy, normal] = ellipsoid3D_equal(scale,level)

% semi-major axis lengths

a = 1.2;
b = 1.2;
c = 0.7;

[~,E] = ellipke(1-(scale*c)^2/(scale*a)^2);
arclength_phi = 2*scale*a*E;

phi = 0;
for i = 1:(level-1)
    F = @(t) (scale*sqrt((a*sin(t)-a*sin(phi(i)))^2+(c*cos(t)-c*cos(phi(i)))^2))-...
        arclength_phi/level;
    t = fzero(F,[phi(i),phi(i)+2*pi/level]);
    phi = [phi,t];
end
phi = [phi,pi];

bdy = [0;0;scale*c];
normal = [0;0;1];
for i = 2:length(phi)-1
    temp = ceil(2*pi*scale*a*sin(phi(i))/arclength_phi*level);

    % this choice of theta only works for a = b

    if mod(i,2) == 1
        theta = 0:2*pi/temp:2*pi-pi/temp;
    else
        theta = (0:2*pi/temp:2*pi-pi/temp)+pi/2;
    end

    % boundary points

    bdy_temp = scale*[a*cos(theta).*sin(phi(i));b*sin(theta).*sin(phi(i));...
        c*cos(phi(i)).*ones(1,length(theta))];
    bdy = cat(2,bdy,bdy_temp);

    % normal vectors

    normalization = sqrt(4./a.^2.*cos(theta).^2.*sin(phi(i)).^2+4./b.^2*sin(theta).^2.*...
        sin(phi(i)).^2+4./c.^2*cos(phi(i)).^2);
    normal_temp = [2./a.*sin(phi(i)).*cos(theta); 2./b.*sin(phi(i)).*sin(theta);...
        2./c.*cos(phi(i)).*ones(1,length(theta))]./...
        normalization;
    normal = cat(2,normal,normal_temp);
end

bdy = cat(2,bdy,[0;0;-scale*0.7]);
normal = cat(2,normal,[0;0;-1]);

end

