clc; clear all; close all;

% Compute KM image for data from a point target. Compute the asymptotic KM
% image for a point target. Then compute the absolute difference between
% the true and asymptotic KM images. All units are with respect to meters
% (m) and seconds (s).
%
% This is the code used to generate the images in the manuscript,
% "Detecting, locating and spectrally characterizing targets with SAR," by
% A. D. Kim, J. Simpson, and C. Tsogka. The parameters in the code below
% are those used to generate Fig. 2 in that manuscript.
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

%% Point target data simulation

% target properties

rho = 1;                        % reflectivity
point_loc = [1;-1;0];           % target location

% SAR data computation

data = zeros(k_num,array_size);

for i = 1:k_num
    for j = 1:array_size
        dist = vecnorm(array(:,j)-point_loc,2);
        data(i,j) = rho/(4*pi)*exp(2i*k(i)*dist)/dist^2;
    end
end

%% Image computations

% True KM image computation

search_x = linspace(-250/k0,250/k0,101);

[X, Y] = meshgrid(search_x);

X = X + point_loc(1);
Y = Y + point_loc(2);

I_KM = 0*X;

for m = 1:length(k)
    for n = 1:array_size
        dist = sqrt((X-array(1,n)).^2+(Y-array(2,n)).^2+(array(3,n)).^2);
        I_KM = I_KM + data(m,n)*exp(-1i.*2.*k(m).*dist);
    end
end

I_KM_norm = abs(I_KM)./max(abs(I_KM(:)));

% True KM image plotting

figure('DefaultAxesFontSize',24)
fig = pcolor(k0*search_x,k0*search_x,I_KM_norm);
set(fig,'EdgeColor','None');
xlabel('$k_0(x-x_0)$', 'Interpreter','Latex');
ylabel('$k_0(y-y_0)$', 'Interpreter','Latex');
colorbar
axis square
hold on
plot(point_loc(1),point_loc(2),'r+','LineWidth',2,'MarkerSize',8);
hold off
% print('-depsc', 'point_target_asymptotic_image1.eps');
% print('-dpdf', 'point_target_asymptotic_image1.pdf');

% Asymptotic KM image computation

t = atan(R/H);
phi0 = (Y-point_loc(2))*sin(t)+(point_loc(1)^2-X.^2)/(2*L)+(point_loc(2)^2-...
    Y.^2)*(1+cos(2*t))/(4*L);

I_KM_asy = rho/(4*pi*L^2)*S_M(2*pi*B*(Y-point_loc(2))/c*R/L,k_num).*...
    S_M(2*pi*a*(X-point_loc(1))/(lambda0*L),array_size).*...
    exp(2i*k0*phi0);

I_KM_asy_norm = abs(I_KM_asy)/max(abs(I_KM_asy(:)));

% Asymptotic KM image plotting

figure('DefaultAxesFontSize',24)
fig = pcolor(k0*search_x,k0*search_x,I_KM_asy_norm);
set(fig,'EdgeColor','None');
xlabel('$k_0(x-x_0)$', 'Interpreter','Latex');
ylabel('$k_0(y-y_0)$', 'Interpreter','Latex');
colorbar
axis square
hold on
plot(point_loc(1),point_loc(2),'r+','LineWidth',2,'MarkerSize',8);
hold off
% print('-depsc', 'point_target_asymptotic_image2.eps');
% print('-dpdf', 'point_target_asymptotic_image2.pdf');

% Asymptotic KM image error plotting

figure('DefaultAxesFontSize',24)
fig = pcolor(k0*search_x,k0*search_x,abs(I_KM_norm-I_KM_asy_norm));
set(fig,'EdgeColor','None');
xlabel('$k_0(x-x_0)$', 'Interpreter','Latex');
ylabel('$k_0(y-y_0)$', 'Interpreter','Latex');
colorbar
axis square
hold on
plot(point_loc(1),point_loc(2),'r+','LineWidth',2,'MarkerSize',8);
hold off
% print('-depsc', 'point_target_asymptotic_image3.eps');
% print('-dpdf', 'point_target_asymptotic_image3.pdf');

% sinc-like function S_n(x)

function S = S_M(x,M)
    array_dim = size(x);
    row_size = array_dim(1);
    col_size = array_dim(2);
    S = zeros(array_dim);
    for i = 1:row_size
        for j = 1:col_size
            if x(i,j) == 0
                S(i,j) = M;
            else
                S(i,j) = sin(M*x(i,j)/(M-1)).*csc(x(i,j)/(M-1));
            end
        end
    end
end
