clc; clear all; close all;

% Compute KM images for data from dielectric spheres. Compute the
% asymptotic KM images for dielectric spheres. Then compute the absolute
% difference between the true and asymptotic KM images. All units are with
% respect to meters (m) and seconds (s).
%
% This is the code used to generate the images in the manuscript,
% "Detecting, locating and spectrally characterizing targets with SAR," by
% A. D. Kim, J. Simpson, and C. Tsogka. The parameters in the code below
% are those used to generate Figs. 3 and 4 in that manuscript.
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

%% Case 1) nrel->inf, ka0->0

% sphere parameters

nrel = 1e16;                            % refractive index of sphere
radius = 0.01*lambda0;                  % radius of sphere

% asymptotic f_B^sca

f_asy = -5/6*k.^2*radius^3;

% MFS for SAR data

delta_s = 0.5*radius;                   
delta_int = 0.5*radius;
points = 100;

[bdy,normal] = sphere3D(radius,points);

data = SAR_data_MFS_3D(bdy,normal,nrel,k,array,delta_s,delta_int);

% True KM image computation

search_x = linspace(-250/k0,250/k0,101);

[X, Y] = meshgrid(search_x);

I_KM = 0*X;

for m = 1:length(k)
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

% True KM image plotting

figure('DefaultAxesFontSize',24)
fig = pcolor(k0*search_x,k0*search_x,I_KM_norm);
set(fig,'EdgeColor','None');
xlabel('$k_0(x-x_0)$', 'Interpreter','Latex');
ylabel('$k_0(y-y_0)$', 'Interpreter','Latex');
colorbar
axis square
hold on
plot(radius*k0*cos(s),radius*k0*sin(s),'w:','LineWidth',2); %boundary
plot(0,0,'rx','LineWidth',1.5,'MarkerSize',8); %center
plot(k0*search_x(y),k0*search_x(x),'r+','LineWidth',1.5,'MarkerSize',8); %KM peak
hold off
% print('-depsc', 'asymptotic_KM_images1.eps');
% print('-dpdf', 'asymptotic_KM_images1.pdf');

% Asymptotic KM image computation

t = atan(R/H);
phi0 = Y*sin(t)-X.^2/(2*L)-Y.^2*(1+cos(2*t))/(4*L);

I_KM_asy = 0*X;

for m = 1:length(k)
    for n = 1:array_size
        tau = Y*sin(t)+array(1,n)*X/L-X.^2/(2*L)-Y.^2*(1+cos(2*t))/(4*L);
        I_KM_asy = I_KM_asy + f_asy(m)*exp(2i*k(m)*tau);
    end
end

I_KM_asy_norm = abs(I_KM_asy)/max(abs(I_KM_asy(:)));

K = abs(I_KM_asy_norm);
[max_num,max_idx] = max(K(:));
[x,y]=ind2sub(size(K),max_idx);

% Asymptotic KM image plotting

figure('DefaultAxesFontSize',24)
fig = pcolor(k0*search_x,k0*search_x,I_KM_asy_norm);
set(fig,'EdgeColor','None');
xlabel('$k_0(x-x_0)$', 'Interpreter','Latex');
ylabel('$k_0(y-y_0)$', 'Interpreter','Latex');
colorbar
axis square
hold on
plot(radius*k0*cos(s),radius*k0*sin(s),'w:','LineWidth',2); %boundary
plot(0,0,'rx','LineWidth',1.5,'MarkerSize',8); %center
plot(k0*search_x(y),k0*search_x(x),'r+','LineWidth',1.5,'MarkerSize',8); %KM peak
hold off
% print('-depsc', 'asymptotic_KM_images2.eps');
% print('-dpdf', 'asymptotic_KM_images2.pdf');

% Asymptotic KM image error plotting

figure('DefaultAxesFontSize',24)
fig = pcolor(k0*search_x,k0*search_x,abs(I_KM_norm-I_KM_asy_norm));
set(fig,'EdgeColor','None');
xlabel('$k_0(x-x_0)$', 'Interpreter','Latex');
ylabel('$k_0(y-y_0)$', 'Interpreter','Latex');
colorbar
axis square
hold on
plot(radius*k0*cos(s),radius*k0*sin(s),'w:','LineWidth',2); %boundary
plot(0,0,'rx','LineWidth',1.5,'MarkerSize',8); %center
plot(k0*search_x(y),k0*search_x(x),'r+','LineWidth',1.5,'MarkerSize',8); %KM peak
hold off
% print('-depsc', 'asymptotic_KM_images3.eps');
% print('-dpdf', 'asymptotic_KM_images3.pdf');

%% Case 2) nrel->inf, ka0->inf

% sphere parameters

nrel = 1e16;                            % refractive index of sphere
radius = 3*lambda0;                     % radius of sphere

% asymptotic f_B^sca

f_asy = -radius/2*exp(-2i*k*radius);

% MFS for SAR data

delta_s = 0.6*radius;
delta_int = 0.6*radius;
points = 800;

[bdy,normal] = sphere3D(radius,points);

data = SAR_data_MFS_3D(bdy,normal,nrel,k,array,delta_s,delta_int);

% True KM image computation

search_x = linspace(-250/k0,250/k0,101);

[X, Y] = meshgrid(search_x);

I_KM = 0*X;

for m = 1:length(k)
    for n = 1:array_size
        dist = sqrt((X-array(1,n)).^2+(Y-array(2,n)).^2+(array(3,n)).^2);
        I_KM = I_KM + data(m,n)*exp(-1i.*2.*k(m).*dist);
    end
end

I_KM_norm = abs(I_KM)./max(abs(I_KM(:)));

K = abs(I_KM_norm);
[max_num,max_idx] = max(K(:));
[x,y]=ind2sub(size(K),max_idx);

% True KM image plotting

s = linspace(0,2*pi,21);
figure('DefaultAxesFontSize',24)
fig = pcolor(k0*search_x,k0*search_x,I_KM_norm);
set(fig,'EdgeColor','None');
xlabel('$k_0(x-x_0)$', 'Interpreter','Latex');
ylabel('$k_0(y-y_0)$', 'Interpreter','Latex');
colorbar
axis square
hold on
plot(radius*k0*cos(s),radius*k0*sin(s),'w:','LineWidth',2); %boundary
plot(0,0,'rx','LineWidth',1.5,'MarkerSize',8); %center
plot(k0*search_x(y),k0*search_x(x),'r+','LineWidth',1.5,'MarkerSize',8); %KM peak
hold off
% print('-depsc', 'asymptotic_KM_images4.eps');
% print('-dpdf', 'asymptotic_KM_images4.pdf');

% Asymptotic KM image computation

t = atan(R/H);
phi0 = Y*sin(t)-X.^2/(2*L)-Y.^2*(1+cos(2*t))/(4*L);

I_KM_asy = radius/2/(4*pi*L^2)*S_M(2*pi*B*(Y-radius*L/R)/c*R/L,k_num).*...
    S_M(2*pi*a*X/(lambda0*L),array_size).*...
    exp(2i*k0*(phi0-radius));

I_KM_asy_norm = abs(I_KM_asy)/max(abs(I_KM_asy(:)));

K = abs(I_KM_asy_norm);
[max_num,max_idx] = max(K(:));
[x,y]=ind2sub(size(K),max_idx);

% Asymptotic KM image plotting

figure('DefaultAxesFontSize',24)
fig = pcolor(k0*search_x,k0*search_x,I_KM_asy_norm);
set(fig,'EdgeColor','None');
xlabel('$k_0(x-x_0)$', 'Interpreter','Latex');
ylabel('$k_0(y-y_0)$', 'Interpreter','Latex');
colorbar
axis square
hold on
plot(radius*k0*cos(s),radius*k0*sin(s),'w:','LineWidth',2); %boundary
plot(0,0,'rx','LineWidth',1.5,'MarkerSize',8); %center
plot(k0*search_x(y),k0*search_x(x),'r+','LineWidth',1.5,'MarkerSize',8); %KM peak
hold off
% print('-depsc', 'asymptotic_KM_images5.eps');
% print('-dpdf', 'asymptotic_KM_images5.pdf');

% Asymptotic KM image error plotting

figure('DefaultAxesFontSize',24)
fig = pcolor(k0*search_x,k0*search_x,abs(I_KM_norm-I_KM_asy_norm));
set(fig,'EdgeColor','None');
xlabel('$k_0(x-x_0)$', 'Interpreter','Latex');
ylabel('$k_0(y-y_0)$', 'Interpreter','Latex');
colorbar
axis square
hold on
plot(radius*k0*cos(s),radius*k0*sin(s),'w:','LineWidth',2); %boundary
plot(0,0,'rx','LineWidth',1.5,'MarkerSize',8); %center
plot(k0*search_x(y),k0*search_x(x),'r+','LineWidth',1.5,'MarkerSize',8); %KM peak
hold off
% print('-depsc', 'asymptotic_KM_images6.eps');
% print('-dpdf', 'asymptotic_KM_images6.pdf');

%% Case 3) nrel->1, ka0->inf

% sphere parameters

nrel = 1.001;                           % refractive index of sphere
radius = 3*lambda0;                     % radius of sphere

% asymptotic f_B^sca

f_asy = -radius/2*exp(-2i*k*radius);

% MFS for SAR data

delta_s = 0.5*radius;
delta_int = 0.5*radius;
points = 1200;

[bdy,normal] = sphere3D(radius,points);

data = SAR_data_MFS_3D(bdy,normal,nrel,k,array,delta_s,delta_int);

% True KM image computation

search_x = linspace(-250/k0,250/k0,101);

[X, Y] = meshgrid(search_x);

I_KM = 0*X;

for m = 1:length(k)
    for n = 1:array_size
        dist = sqrt((X-array(1,n)).^2+(Y-array(2,n)).^2+(array(3,n)).^2);
        I_KM = I_KM + data(m,n)*exp(-1i.*2.*k(m).*dist);
    end
end

I_KM_norm = abs(I_KM)./max(abs(I_KM(:)));

K = abs(I_KM_norm);
[max_num,max_idx] = max(K(:));
[x,y]=ind2sub(size(K),max_idx);

% True KM image plotting

s = linspace(0,2*pi,21);
figure('DefaultAxesFontSize',24)
fig = pcolor(k0*search_x,k0*search_x,I_KM_norm);
set(fig,'EdgeColor','None');
xlabel('$k_0(x-x_0)$', 'Interpreter','Latex');
ylabel('$k_0(y-y_0)$', 'Interpreter','Latex');
colorbar
axis square
hold on
plot(radius*k0*cos(s),radius*k0*sin(s),'w:','LineWidth',2); %boundary
plot(0,0,'rx','LineWidth',1.5,'MarkerSize',8); %center
plot(k0*search_x(y),k0*search_x(x),'r+','LineWidth',1.5,'MarkerSize',8); %KM peak
hold off
% print('-depsc', 'asymptotic_KM_images7.eps');
% print('-dpdf', 'asymptotic_KM_images7.pdf');

% Asymptotic KM Image

t = atan(R/H);
phi0 = Y*sin(t)-X.^2/(2*L)-Y.^2*(1+cos(2*t))/(4*L);

I_KM_asy = (radius*(nrel-1)/(4*nrel))/(4*pi*L^2).*(S_M(2*pi*B*(Y+radius*(nrel+1)*L/(2*R))/c*R/L,k_num).*...
    S_M(2*pi*a*X/(lambda0*L),array_size).*exp(2i*k0*(phi0+radius*(nrel+1)/2))+...
    S_M(2*pi*B*(Y-radius*(nrel+1)*L/(2*R))/c*R/L,k_num).*...
    S_M(2*pi*a*X/(lambda0*L),array_size).*exp(2i*k0*(phi0-radius*(nrel+1)/2)));

I_KM_asy_norm = abs(I_KM_asy)/max(abs(I_KM_asy(:)));

K = abs(I_KM_asy_norm);
[max_num,max_idx] = max(K(:));
[x,y]=ind2sub(size(K),max_idx);

% Asymptotic KM image computation

figure('DefaultAxesFontSize',24)
fig = pcolor(k0*search_x,k0*search_x,I_KM_asy_norm);
set(fig,'EdgeColor','None');
xlabel('$k_0(x-x_0)$', 'Interpreter','Latex');
ylabel('$k_0(y-y_0)$', 'Interpreter','Latex');
colorbar
axis square
hold on
plot(radius*k0*cos(s),radius*k0*sin(s),'w:','LineWidth',2); %boundary
plot(0,0,'rx','LineWidth',1.5,'MarkerSize',8); %center
plot(k0*search_x(y),k0*search_x(x),'r+','LineWidth',1.5,'MarkerSize',8); %KM peak
hold off
% print('-depsc', 'asymptotic_KM_images8.eps');
% print('-dpdf', 'asymptotic_KM_images8.pdf');

% Asymptotic KM image error plotting

figure('DefaultAxesFontSize',24)
fig = pcolor(k0*search_x,k0*search_x,abs(I_KM_norm-I_KM_asy_norm));
set(fig,'EdgeColor','None');
xlabel('$k_0(x-x_0)$', 'Interpreter','Latex');
ylabel('$k_0(y-y_0)$', 'Interpreter','Latex');
colorbar
axis square
hold on
plot(radius*k0*cos(s),radius*k0*sin(s),'w:','LineWidth',2); %boundary
plot(0,0,'rx','LineWidth',1.5,'MarkerSize',8); %center
plot(k0*search_x(y),k0*search_x(x),'r+','LineWidth',1.5,'MarkerSize',8); %KM peak
hold off
% print('-depsc', 'asymptotic_KM_images9.eps');
% print('-dpdf', 'asymptotic_KM_images9.pdf');

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