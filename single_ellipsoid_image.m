clc; clear all; close all;

% Compute KM image for data from an ellipsoid target. Then compute f^sca_B.
% All units are with respect to meters (m) and seconds (s).
%
% This is the code used to generate the images in the manuscript,
% "Spectrally characterizing targets with SAR," by
% A. D. Kim, J. Simpson, and C. Tsogka. The parameters in the code below
% are those used to generate Fig. 5 in that manuscript.
%
% Written by J. Simpson on 5/8/2025

%% Imaging System Parameters

% bandwidth

k_num = 26;                     % number of frequency samples

B = 622e6;                      % system bandwidth (Hz)
c = 3e8;                        % wave speed (m/s)
f0 = 9.6e9;                     % central frequency (Hz)
k0 = 2*pi*f0/c;                 % central wavenumber (m^-1)
dk = pi*B/c;                    % maximum wavenumber distance (m^-1)

k = linspace(k0-dk,k0+dk,k_num);% wavenumber array

lambda0 = 0.0312;               % central wavelength

% flight path

L = 8120;                       % distance from array to imaging window (m)
R = 3550;                       % range offset (m)
H = 7300;                       % height (m)
a = 130;                        % aperture length (m)

array_size = 32;                % number of spatial samples

array = zeros(3,array_size);    % flight path array

for j = 1:array_size
    array(1,j)=(a*(2*j-array_size-1))/(2*array_size-2);
end
array(2,:) = ones(1,array_size)*(R);
array(3,:) = ones(1,array_size)*(H);

r = pi/4;
rot = [[cos(r),0,sin(r)];
    [0,1,0];
    [-sin(r),0,cos(r)]];

array = rot*array;

%% SAR data simulation

% target parameters

nrel = 1e16;                        % refractive index of ellipsoid
scale = 2.5*lambda0;                % size of ellipsoid

% MFS for SAR data

delta_s = 0.3*scale;
delta_int = 0.3*scale;
level = 25;

[bdy,normal] = ellipsoid3D_equal(scale,level);

data = SAR_data_MFS_3D(bdy,normal,nrel,k,array,delta_s,delta_int);

%% KM image computation

search_x = linspace(-250/k0,250/k0,501);

array = inv(rot)*array;

[X, Y] = meshgrid(search_x);

I_KM = 0*X;

for m = 1:k_num
    for n = 1:array_size
        dist = sqrt((X-array(1,n)).^2+(Y-array(2,n)).^2+(array(3,n)).^2);
        I_KM = I_KM + data(m,n)*exp(-1i.*2.*k(m).*dist);
    end
end

I_KM_norm = abs(I_KM)./max(abs(I_KM(:)));

K = abs(I_KM_norm);
[max_num,max_idx] = max(K(:));
[x,y]=ind2sub(size(K),max_idx);

s = linspace(0,2*pi,21);

% KM image plotting

figure('DefaultAxesFontSize',24)
fig = pcolor(k0*search_x,k0*search_x,I_KM_norm);
set(fig,'EdgeColor','None');
xlabel('$k_0(x-x_0)$', 'Interpreter','Latex');
ylabel('$k_0(y-y_0)$', 'Interpreter','Latex');
colorbar
axis square
hold on
plot(k0*search_x(y),k0*search_x(x),'r+','LineWidth',1.5); %peak location
plot(0,0,'rx','LineWidth',1.5);
hold off
% print('-depsc', 'single_ellipsoid_image1.eps');
% print('-dpdf', 'single_ellipsoid_image1.pdf');

%% Scattering amplitude computation

dir = rot*[0;R;H];

dir = dir/vecnorm(dir,2);

f1 = backscatamp3D_mfs(bdy,normal,delta_s,delta_int,nrel,k,dir);

f2 = zeros(1,21);
for i = 1:21
    dir1 = rot*[-a/2+(i-1)/20*a;R;H];
    dir1 = dir1/vecnorm(dir1,2);
    f2(i) = backscatamp3D_mfs(bdy,normal,delta_s,delta_int,nrel,k0,dir1);
end

% Scattering amplitude plotting

freq = c*k/(2*pi)*1e-9;
phi1 = 2*atan(a/(2*L))*180/pi;
phi = linspace(-phi1,phi1,21);

figure('DefaultAxesFontSize',20)
plot(freq,real(f1),'k','LineWidth',2); hold on
plot(freq,imag(f1),'r','LineWidth',2);
xlim([min(freq),max(freq)]);
xlabel('frequency (GHz)', 'Interpreter','Latex');
ylabel('$f^{sca}_B$ (m)', 'Interpreter','Latex');
legend('real','imaginary',...
    'Interpreter','Latex','Location','Best');
% print('-depsc', 'single_ellipsoid_image2.eps');
% print('-dpdf', 'single_ellipsoid_image2.pdf');

figure('DefaultAxesFontSize',20)
plot(phi,real(f2),'k','LineWidth',2); hold on
plot(phi,imag(f2),'r','LineWidth',2);
xlim([min(phi),max(phi)]);
xlabel('$\phi$ (deg)', 'Interpreter','Latex');
ylabel('$f^{sca}_B$ (m)', 'Interpreter','Latex');
legend('real','imaginary',...
    'Interpreter','Latex','Location','Best');
% print('-depsc', 'single_ellipsoid_image3.eps');
% print('-dpdf', 'single_ellipsoid_image3.pdf');
