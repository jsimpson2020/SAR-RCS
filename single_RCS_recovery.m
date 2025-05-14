clc; clear all; close all;

% Compute KM images for data from dielectric finite-sized targets. Then
% compute the recovered RCS spectra. All units are with respect to meters
% (m) and seconds (s).
%
% This is the code used to generate the images in the manuscript,
% "Spectrally characterizing targets with SAR," by
% A. D. Kim, J. Simpson, and C. Tsogka. The parameters in the code below
% are those used to generate Fig. 7.
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

%% SAR data simulation

% target parameters

nrel = 1.4;                         % refractive index
radius = 2.5*lambda0;               % radius

% MFS for SAR data

delta_s = 0.4*radius;
delta_int = 0.4*radius;
points = 800;

[bdy,normal] = sphere3D(radius,points);

data = SAR_data_MFS_3D(bdy,normal,nrel,k,array,delta_s,delta_int);

%% Add measurement noise
  
sn     = 0.1; % this is the parameter to change
B0     = data;
enr    = sqrt(sum(abs(data(:)).^2)/length(data(:)));
Bnoise = enr*sn*(randn(size(data))+1i*randn(size(data)))/sqrt(2);
BB     = B0 + Bnoise;

SNR_est       = -10 * log10( sn );
SNR_real      =  10 * log10( norm(data) / norm(Bnoise) );
   
disp( ' ' );
disp( '   ---' );
disp( '   SNR' );
disp( '   ---' );
disp( ' ' );
disp( ['   Realization SNR = ', num2str(SNR_real), ' dB' ] );
disp( ['   Average SNR     = ', num2str(SNR_est), ' dB'  ] );
disp( ' ' );

%% KM image computation

search_x = linspace(-250/k0,250/k0,101);

[X, Y] = meshgrid(search_x);

I_KM = 0*X;

for m = 1:length(k)
    for n = 1:array_size
        dist = sqrt((X-array(1,n)).^2+(Y-array(2,n)).^2+(array(3,n)).^2);
        I_KM = I_KM + BB(m,n)*exp(-1i.*2.*k(m).*dist);
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
plot(radius*k0*cos(s),radius*k0*sin(s),'w:','LineWidth',2); %boundary
plot(0,0,'rx','LineWidth',1.5,'MarkerSize',8); %center
plot(k0*search_x(y),k0*search_x(x),'r+','LineWidth',1.5,'MarkerSize',8); %KM peak
hold off
% print('-depsc', 'single_RCS_recovery1.eps');
% print('-dpdf', 'single_RCS_recovery1.pdf');

%% RCS Recovery

ystar = [search_x(y);search_x(x);0];

F = 0*BB;

for i = 1:array_size
    F(:,i) = vecnorm(array(:,i)-ystar,2)^4.*abs(BB(:,i)).^2;
end
F = (4*pi)^2*F;

sigma_RCS = 4*pi/array_size*sum(F,2);

trunc = 60;

sigma_true = 4*pi*abs(backscatamp3D_anl(trunc,nrel,k,radius)).^2;

freq = c*k/(2*pi)*1e-9;

% RCS plotting

figure('DefaultAxesFontSize',20)
plot(freq,sigma_true,'k','LineWidth',2); hold on
plot(freq,sigma_RCS,'ro','LineWidth',2);
xlim([min(freq),max(freq)]);
xlabel('frequency (GHz)', 'Interpreter','Latex');
ylabel('$\sigma_{RCS}$ (m$^2$)', 'Interpreter','Latex');
legend('true','recovered',...
    'Interpreter','Latex','Location','Best');
% print('-depsc', 'single_RCS_recovery2.eps');
% print('-dpdf', 'single_RCS_recovery2.pdf');

%% SAR data simulation

% target parameters

nrel = 1.6;                         % refractive index
radius = 2.5*lambda0;               % radius

% MFS for SAR data

delta_s = 0.4*radius;
delta_int = 0.4*radius;
points = 1000;

[bdy,normal] = sphere3D(radius,points);

data = SAR_data_MFS_3D(bdy,normal,nrel,k,array,delta_s,delta_int);

%% Add measurement noise
  
sn     = 0.1; % this is the parameter to change
B0     = data;
enr    = sqrt(sum(abs(data(:)).^2)/length(data(:)));
Bnoise = enr*sn*(randn(size(data))+1i*randn(size(data)))/sqrt(2);
BB     = B0 + Bnoise;

SNR_est       = -10 * log10( sn );
SNR_real      =  10 * log10( norm(data) / norm(Bnoise) );
   
disp( ' ' );
disp( '   ---' );
disp( '   SNR' );
disp( '   ---' );
disp( ' ' );
disp( ['   Realization SNR = ', num2str(SNR_real), ' dB' ] );
disp( ['   Average SNR     = ', num2str(SNR_est), ' dB'  ] );
disp( ' ' );

%% KM image computation

search_x = linspace(-250/k0,250/k0,101);

[X, Y] = meshgrid(search_x);

I_KM = 0*X;

for m = 1:length(k)
    for n = 1:array_size
        dist = sqrt((X-array(1,n)).^2+(Y-array(2,n)).^2+(array(3,n)).^2);
        I_KM = I_KM + BB(m,n)*exp(-1i.*2.*k(m).*dist);
    end
end

I_KM_norm = abs(I_KM)./max(abs(I_KM(:)));

K = abs(I_KM_norm);
[max_num,max_idx] = max(K(:));
[x,y]=ind2sub(size(K),max_idx);

% KM image plotting

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
% print('-depsc', 'single_RCS_recovery3.eps');
% print('-dpdf', 'single_RCS_recovery3.pdf');

%% RCS Recovery

ystar = [search_x(y);search_x(x);0];

F = 0*BB;

for i = 1:array_size
    F(:,i) = vecnorm(array(:,i)-ystar,2)^4.*abs(BB(:,i)).^2;
end
F = (4*pi)^2*F;

sigma_RCS = 4*pi/array_size*sum(F,2);

trunc = 60;

sigma_true = 4*pi*abs(backscatamp3D_anl(trunc,nrel,k,radius)).^2;

freq = c*k/(2*pi)*1e-9;

% RCS plotting

figure('DefaultAxesFontSize',20)
plot(freq,sigma_true,'k','LineWidth',2); hold on
plot(freq,sigma_RCS,'ro','LineWidth',2);
xlim([min(freq),max(freq)]);
xlabel('frequency (GHz)', 'Interpreter','Latex');
ylabel('$\sigma_{RCS}$ (m$^2$)', 'Interpreter','Latex');
legend('true','recovered',...
    'Interpreter','Latex','Location','Best');
% print('-depsc', 'single_RCS_recovery4.eps');
% print('-dpdf', 'single_RCS_recovery4.pdf');

%% SAR data simulation

% target parameters

nrel = 1.4;                         % refractive index
scale = 2.5*lambda0;                % radius

% MFS for SAR data

delta_s = 0.3*scale;
delta_int = 0.3*scale;
level = 25;

[bdy,normal] = ellipsoid3D_equal(scale,level);

data = SAR_data_MFS_3D(bdy,normal,nrel,k,array,delta_s,delta_int);

%% Add measurement noise
  
sn     = 0.1; % this is the parameter to change
B0     = data;
enr    = sqrt(sum(abs(data(:)).^2)/length(data(:)));
Bnoise = enr*sn*(randn(size(data))+1i*randn(size(data)))/sqrt(2);
BB     = B0 + Bnoise;

SNR_est       = -10 * log10( sn );
SNR_real      =  10 * log10( norm(data) / norm(Bnoise) );
   
disp( ' ' );
disp( '   ---' );
disp( '   SNR' );
disp( '   ---' );
disp( ' ' );
disp( ['   Realization SNR = ', num2str(SNR_real), ' dB' ] );
disp( ['   Average SNR     = ', num2str(SNR_est), ' dB'  ] );
disp( ' ' );

%% KM image computation

search_x = linspace(-250/k0,250/k0,101);

[X, Y] = meshgrid(search_x);

I_KM = 0*X;

for m = 1:length(k)
    for n = 1:array_size
        dist = sqrt((X-array(1,n)).^2+(Y-array(2,n)).^2+(array(3,n)).^2);
        I_KM = I_KM + BB(m,n)*exp(-1i.*2.*k(m).*dist);
    end
end

I_KM_norm = abs(I_KM)./max(abs(I_KM(:)));

K = abs(I_KM_norm);
[max_num,max_idx] = max(K(:));
[x,y]=ind2sub(size(K),max_idx);

% KM image plotting

figure('DefaultAxesFontSize',24)
fig = pcolor(k0*search_x,k0*search_x,I_KM_norm);
set(fig,'EdgeColor','None');
xlabel('$k_0(x-x_0)$', 'Interpreter','Latex');
ylabel('$k_0(y-y_0)$', 'Interpreter','Latex');
colorbar
axis square
hold on
plot(1.2*radius*k0*cos(s),1.2*radius*k0*sin(s),'w:','LineWidth',2); %boundary
plot(0,0,'rx','LineWidth',1.5,'MarkerSize',8); %center
plot(k0*search_x(y),k0*search_x(x),'r+','LineWidth',1.5,'MarkerSize',8); %KM peak
hold off
% print('-depsc', 'single_RCS_recovery5.eps');
% print('-dpdf', 'single_RCS_recovery5.pdf');

%% RCS Recovery

ystar = [search_x(y);search_x(x);0];

F = 0*BB;

for i = 1:array_size
    F(:,i) = vecnorm(array(:,i)-ystar,2)^4.*abs(BB(:,i)).^2;
end
F = (4*pi)^2*F;

sigma_RCS = 4*pi/array_size*sum(F,2);

dir = [0;R;H];
dir = dir/vecnorm(dir,2);

sigma_true = 4*pi*abs(backscatamp3D_mfs(bdy,normal,delta_s,delta_int,nrel,k,dir)).^2;

freq = c*k/(2*pi)*1e-9;

% RCS plotting

figure('DefaultAxesFontSize',20)
plot(freq,sigma_true,'k','LineWidth',2); hold on
plot(freq,sigma_RCS,'ro','LineWidth',2);
xlim([min(freq),max(freq)]);
xlabel('frequency (GHz)', 'Interpreter','Latex');
ylabel('$\sigma_{RCS}$ (m$^2$)', 'Interpreter','Latex');
legend('true','recovered',...
    'Interpreter','Latex','Location','Best');
% print('-depsc', 'single_RCS_recovery6.eps');
% print('-dpdf', 'single_RCS_recovery6.pdf');
