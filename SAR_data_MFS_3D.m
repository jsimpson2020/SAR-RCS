% This function computes SAR data given boundary points and normal vectors
% for a penetrable target using the Method of Fundamental Solutions (MFS)
%
% bdy = boundary points of target
% normal = normal vectors for boundary points of target
% nrel = relative refractive index of target
% k = wavenumber vector
% array = imaging array vector
% delta_s = scattered field delta for MFS
% delta_int = interior field delta for MFS
%
% Written by J. Simpson on 5/8/2025

function data = SAR_data_MFS_3D(bdy,normal,nrel,k,array,delta_s,delta_int)

M = length(k);                  % data matrix row size
N = length(array(1,:));         % data matrix col size
points = length(bdy(1,:));      % number of MFS points

data = zeros(M,N);              % data matrix

source_scat = bdy - delta_s.* normal; % MFS source points for scattered field
source_int = bdy + delta_int.* normal;% MFS source points for interior field

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

rhs_dist = zeros(points,N);
rhs_normdot = zeros(points,N);
sar_dist = zeros(points,N);
for i = 1:points
    for j = 1:N
        rhs_dist(i,j) = sqrt((bdy(:,i)-array(:,j))'*(bdy(:,i)-array(:,j)));
        rhs_normdot(i,j) = normal(:,i)'*(bdy(:,i)-array(:,j));
        sar_dist(i,j) = sqrt((source_scat(:,i)-array(:,j))'*(source_scat(:,i)-array(:,j)));
    end
end

% Build and solve MFS linear system to compute data

for i = 1:M
    A = exp(1i*k(i)*dist_scat)./(4*pi*dist_scat);
    B = -exp(1i*nrel.*k(i).*dist_int)./(4*pi*dist_int);
    C = (-1+1i*k(i)*dist_scat).*exp(1i*k(i)*dist_scat).*normdot_scat./(4*pi*dist_scat.^3);
    D = (1./nrel.^2-1i*k(i)./nrel.*dist_int).*exp(1i.*nrel*k(i).*dist_int).*normdot_int./(4*pi*dist_int.^3);

    MFS_matcell = {A, B; C, D};
    MFS_matrix = cell2mat(MFS_matcell);
    
    for j = 1:N
        rhs_1 = -exp(1i*k(i)*rhs_dist(:,j))./(4*pi*rhs_dist(:,j));
        rhs_2 = (1-1i*k(i)*rhs_dist(:,j)).*exp(1i*k(i)*rhs_dist(:,j)).*rhs_normdot(:,j)./(4*pi*rhs_dist(:,j).^3);

        MFS_rhs = [rhs_1; rhs_2];

        c = MFS_matrix\MFS_rhs;

        MFS_weights = c(1:points);

        data(i,j) = MFS_weights.'*(exp(1i*k(i)*sar_dist(:,j))./(4*pi*sar_dist(:,j)));
    end
end

end



