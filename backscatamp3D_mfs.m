% This function computes the scattering amplitude in the direction of
% backscattering given boundary points and normal vectors using the Method
% of Fundamental Solutions (MFS)
%
% bdy = boundary points of target
% normal = normal vectors for boundary points of target
% delta_s = scattered field delta for MFS
% delta_int = interior field delta for MFS
% nrel = relative refractive index of target
% k = wavenumber vector
% dir = incident plane wave direction vector
%
% Written by J. Simpson on 5/8/2025

function fsca = backscatamp3D_mfs(bdy,normal,delta_s,delta_int,nrel,k,dir)

points = length(bdy(1,:));              % number of MFS points

fsca = zeros(length(k),length(nrel));   % scattering amplitude array

source_scat = bdy - delta_s * normal;   % MFS source points for scattered field
source_int = bdy + delta_int * normal;  % MFS source points for interior field

% MFS distance computations

dist_scat = zeros(points,points);
dist_int = zeros(points,points);
normdot_scat = zeros(points,points);
normdot_int = zeros(points,points);
for i = 1:points
    for j = 1:points
        dist_scat(i,j) = sqrt((bdy(:,i)-source_scat(:,j))'*(bdy(:,i)-source_scat(:,j)));
        dist_int(i,j) = sqrt((bdy(:,i)-source_int(:,j))'*(bdy(:,i)-source_int(:,j)));
        normdot_scat(i,j) = normal(:,i)'*(bdy(:,i)-source_scat(:,j));
        normdot_int(i,j) = normal(:,i)'*(bdy(:,i)-source_int(:,j));
    end
end

x_dot_d = dir.'*bdy;
n_dot_d = dir.'*normal;
source_dot_d = dir.'*source_scat;

for i = 1:length(k)
    rhs_1 = exp(1i*k(i)*x_dot_d).';
    rhs_2 = (n_dot_d*1i*k(i).*exp(1i*k(i)*x_dot_d)).';

    MFS_rhs = [rhs_1; rhs_2];

    % Build and solve MFS linear system
    
    for j = 1:length(nrel)
        A = exp(1i*k(i)*dist_scat)./(4*pi*dist_scat);
        B = -exp(1i*nrel(j)*k(i)*dist_int)./(4*pi*dist_int);
        C = (-1+1i*k(i)*dist_scat).*exp(1i*k(i)*dist_scat).*normdot_scat./(4*pi*dist_scat.^3);
        D = (1/nrel(j)^2-1i*k(i)/nrel(j)*dist_int).*exp(1i*nrel(j)*k(i)*dist_int).*normdot_int./(4*pi*dist_int.^3);

        MFS_matcell = {A, B; C, D};
        MFS_matrix = cell2mat(MFS_matcell);

        c = MFS_matrix\MFS_rhs;

        MFS_weights = c(1:points);

        fsca(i,j) = 0.25/pi*exp(1i*k(i)*source_dot_d)*MFS_weights;

    end
end