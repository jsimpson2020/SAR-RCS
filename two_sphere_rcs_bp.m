clc; clear all; close all;

% Compute KM images for data from two dielectric spheres. Then compute the
% recovered RCS spectra. All units are with respect to meters (m) and
% seconds (s).
%
% This is the code used to generate the images in the manuscript,
% "Spectrally characterizing targets with SAR," by
% A. D. Kim, J. Simpson, and C. Tsogka. The parameters in the code below
% are those used to generate Fig. 8.
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

nrel1 = 1.4;                            % refractive index of sphere 1
radius1 = 2*lambda0;                    % radius of sphere 1
center1 = [-40*lambda0;2*lambda0;0];    % center of sphere 1

nrel2 = 1.6;                            % refractive index of sphere 2
radius2 = 1.8*lambda0;                  % radius of sphere 2
center2 = [50*lambda0;-2*lambda0;0];    % center of sphere 2

% MFS parameters for sphere 1

delta_s1 = 0.4*radius1;
delta_int1 = 0.4*radius1;
points1 = 500;

% MFS parameters for sphere 2

delta_s2 = 0.4*radius2;
delta_int2 = 0.4*radius2;
points2 = 500;

% MFS for SAR data

[bdy,normal] = twospheres(radius1,radius2,center1,center2,points1,points2);

nrel1_array = nrel1*ones(1,points1);
nrel2_array = nrel2*ones(1,points2);
nrel = [nrel1_array,nrel2_array];

delta_s1_array = delta_s1*ones(1,points1);
delta_s2_array = delta_s2*ones(1,points2);
delta_s = [delta_s1_array,delta_s2_array];

delta_int1_array = delta_int1*ones(1,points1);
delta_int2_array = delta_int2*ones(1,points2);
delta_int = [delta_int1_array,delta_int2_array];

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

search_x = linspace(-500/k0,500/k0,1001);

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
[x_KM,y_KM]=ind2sub(size(K),max_idx);

s = linspace(0,2*pi,21);

peak1 = [search_x(y_KM);search_x(x_KM);0];

disp(peak1/lambda0)

% KM image plotting

figure('DefaultAxesFontSize',24)
fig = pcolor(k0*search_x,k0*search_x,I_KM_norm);
set(fig,'EdgeColor','None');
xlabel('$k_0x$', 'Interpreter','Latex');
ylabel('$k_0y$', 'Interpreter','Latex');
colorbar
axis square
hold on
plot(k0*search_x(y_KM),k0*search_x(x_KM),'r+','LineWidth',1.5); %peak location
plot(k0*center1(1),k0*center1(2),'rx','LineWidth',1.5);
plot(k0*center1(1)+k0*radius1*cos(s),k0*center1(2)+k0*radius1*sin(s),'w:',...
    'LineWidth',2);
plot(k0*center2(1),k0*center2(2),'rx','LineWidth',1.5);
plot(k0*center2(1)+k0*radius2*cos(s),k0*center2(2)+k0*radius2*sin(s),'w:',...
    'LineWidth',2);
hold off

%% KM Linearity Subtraction Method

Dprime = zeros(k_num,array_size);
tar_loc = [search_x(y_KM);search_x(x_KM);0];

for i = 1:k_num
    for j = 1:array_size
        Dprime(i,j) = exp(2*1i*k(i)*vecnorm(array(:,j)-tar_loc,2))/...
            (vecnorm(array(:,j)-tar_loc,2))^2;
    end
end

I_KM_prime = 0*X;

for m = 1:length(k)
    for n = 1:array_size
        R = sqrt((X-array(1,n)).^2+(Y-array(2,n)).^2+(array(3,n)).^2);
        I_KM_prime = I_KM_prime + Dprime(m,n)*exp(-1i.*2.*k(m).*R);
    end
end

I_KM_prime = abs(I_KM_prime)./max(abs(I_KM_prime(:)));

I_KM_2 = I_KM_norm - I_KM_prime;

I_KM_2_norm = I_KM_2./max(abs(I_KM_2(:)));


% Modified KM image plotting

K = I_KM_2_norm;
[max_num,max_idx] = max(K(:));
[x,y]=ind2sub(size(K),max_idx);

peak2 = [search_x(y);search_x(x);0];

disp(peak2/lambda0)

figure('DefaultAxesFontSize',24)
fig = pcolor(k0*search_x,k0*search_x,I_KM_2_norm);
set(fig,'EdgeColor','None');
xlabel('$k_0x$', 'Interpreter','Latex');
ylabel('$k_0y$', 'Interpreter','Latex');
colorbar
axis square
hold on
plot(k0*search_x(y),k0*search_x(x),'r+','LineWidth',1.5); %peak location
plot(k0*center1(1),k0*center1(2),'rx','LineWidth',1.5);
plot(k0*center1(1)+k0*radius1*cos(s),k0*center1(2)+k0*radius1*sin(s),'w:',...
    'LineWidth',2);
plot(k0*center2(1),k0*center2(2),'rx','LineWidth',1.5);
plot(k0*center2(1)+k0*radius2*cos(s),k0*center2(2)+k0*radius2*sin(s),'w:',...
    'LineWidth',2);
hold off


%% Array Back-propagation Method

f1 = zeros(k_num,1);
f2 = zeros(k_num,1);

diag1 = zeros(k_num,1);
diag2 = zeros(k_num,1);

for i = 1:k_num
    f1_temp = 0;
    f2_temp = 0;
    d1_temp = 0;
    d2_temp = 0;
    for j = 1:array_size
        f1_temp = f1_temp + BB(i,j)*vecnorm(array(:,j)-peak1,2)^2*...
            exp(-2i*k(i)*vecnorm(array(:,j)-peak1,2));
        f2_temp = f2_temp + BB(i,j)*vecnorm(array(:,j)-peak2,2)^2*...
            exp(-2i*k(i)*vecnorm(array(:,j)-peak2,2));
        d1_temp = d1_temp + vecnorm(array(:,j)-peak1,2)^2/...
            vecnorm(array(:,j)-peak2,2)^2*...
            exp(2i*k(i)*(vecnorm(array(:,j)-peak2,2)-vecnorm(array(:,j)-peak1,2)));
        d2_temp = d2_temp + vecnorm(array(:,j)-peak2,2)^2/...
            vecnorm(array(:,j)-peak1,2)^2*...
            exp(2i*k(i)*(vecnorm(array(:,j)-peak1,2)-vecnorm(array(:,j)-peak2,2)));
    end
    f1(i) = f1_temp;
    f2(i) = f2_temp;
    diag1(i) = d1_temp;
    diag2(i) = d2_temp;
end
f1 = 4*pi/array_size*f1;
f2 = 4*pi/array_size*f2;
diag1 = diag1/array_size;
diag2 = diag2/array_size;

f = [f1;f2];

A = [[eye(k_num),diag(diag1)];
    [diag(diag2),eye(k_num)]];

fvec = A\f;

RCS1_rec = 4*pi*abs(fvec(1:k_num)).^2;
RCS2_rec = 4*pi*abs(fvec(k_num+1:end)).^2;

trunc = 60;

RCS1_true = 4*pi*abs(backscatamp3D_anl(trunc,nrel1,k,radius1)).^2;
RCS2_true = 4*pi*abs(backscatamp3D_anl(trunc,nrel2,k,radius2)).^2;

% RCS plotting

freq = c*k/(2*pi)*1e-9;

figure('DefaultAxesFontSize',20)
plot(freq,RCS1_true,'k','LineWidth',2); hold on
plot(freq,RCS1_rec,'ro','LineWidth',2);
xlim([min(freq),max(freq)]);
xlabel('frequency (GHz)', 'Interpreter','Latex');
ylabel('$\sigma_{RCS}$ (m$^2$)', 'Interpreter','Latex');
legend('true','recovered',...
    'Interpreter','Latex','Location','Best');
% print('-depsc', 'two_sphere_rcs_bp1.eps');
% print('-dpdf', 'two_sphere_rcs_bp1.pdf');

figure('DefaultAxesFontSize',20)
plot(freq,RCS2_true,'k','LineWidth',2); hold on
plot(freq,RCS2_rec,'ro','LineWidth',2);
xlim([min(freq),max(freq)]);
xlabel('frequency (GHz)', 'Interpreter','Latex');
ylabel('$\sigma_{RCS}$ (m$^2$)', 'Interpreter','Latex');
legend('true','recovered',...
    'Interpreter','Latex','Location','Best');
% print('-depsc', 'two_sphere_rcs_bp2.eps');
% print('-dpdf', 'two_sphere_rcs_bp2.pdf');

