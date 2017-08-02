%Script to compare corelations between w'T' planer fit from the tower to
%the sensible heat flrux calculated from the FLIR
%
%By: Travis Morrison
%edited 7-28-17
clear; clc; close all;
%% OPTIONS
data_load_option = 'TIR'; %TIR_and_Tower  
%...Load TIR_and_Tower you want comparison with EC method from tower
% Load just TIR if you want just SEB-TIR method 
%data paths
TIR_data_path = '/TIR_data/IOP9/surf_temp1.mat';
Radiation_data_path = 'Radiation_data/playa_05_2013_5minDATA.mat';
%Sonic_data_path = '/PlayaSpringRawData24_28.mat';
%Ground_HF_data_path = '/Users/tjmorrison/Desktop/Reseach/Matterhorn/groundHeatFluxdata/GHF_playa2013spring.mat';

chunk_15 = 1; %decides which 15 minute chunk ~0 correlates to first chunk
%Averaging Window-smooths signal ~reduces noise

Avg_chunk = 10; %gets redefined as tau soil

machine = 'imac_local'; %data path for which machine

%% Preparing data paths and setting plot settings
%add paths
switch machine
    case 'imac_local'
        %code
        addpath('/Users/travismorrison/Documents/Code/functions_library');
        %data path
        data_path ='/Users/travismorrison/Local_Data/MATERHORN/data/';
        addpath(data_path); 
    case 'imac_Achilles'
        %code
        addpath('/Users/travismorrison/Documents/Code/functions_library');
        %data path
        data_path ='/Volumes/ACHILLES/MATERHORN/data/';
        addpath(data_path); 
    case 'linux'
        %code
        addpath('/home/tjmorrison/Research/function_library');
        %data path
        data_path = addpath('/media/tjmorrison/ACHILLES/MATERHORN/data');
end
% Plot xire and TIR Spectras
ft_size = 25;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex'); 
set(0,'DefaultAxesFontSize',ft_size); 
%% Preparing Data

%Ground heat flux data measured
%IOP9~q_groundData.GHF(7164:7176);
%IOP10~ q_groundData.GHF(8644:8684);
%q_groundData = load(Ground_HF_data_path); 
%q_soilData = q_groundData.GHF(7164+chunk_15*3:7176); 


%Grabbing right set of data
%Radition
%IOP 10 [W m^-2]  row 8643, column 10-->NetRad = May 31th 00:15:55.55 UTC
%IOP9 [W m^-2]  row 7164:7176, column 10-->NetRad = May 25th 21:00-22:00 UTC
load(strcat(data_path,Radiation_data_path)) ;
Net_Radiation = playa0520135min(7164+chunk_15*3:7179,10);
SW_in = playa0520135min(7164+chunk_15*3:7179,6);
SW_out = playa0520135min(7164+chunk_15*3:7179,7);
LW_in = playa0520135min(7164+chunk_15*3:7179,8);
LW_out = playa0520135min(7164+chunk_15*3:7179,9);
clear playa0520135min;

switch data_load_option
    case 'TIR'
        %Thermal Camera 
        load(strcat(data_path,TIR_data_path));
        T_surface = surf_temp1;

    case 'Tower_and_TIR'
        
        %Thermal Camera
        load(strcat(data_path,TIR_data_path));
        T_surface = surf_temp1;
        %Tower Data -----> want to revamp with functions
        %IOP9 3236858 correlates to time of May 25th 2013 20:57:22.9 Camera started at 20:57:22.92 UTC
        %IOP10 1747111 correlates to time of May 31th 2013 00:15:55.55 Camera started at 00:15:55.537 UTC
        
        %Sonics
        Sonic_data=load(strcat(data_path,'tower_data/Playa_tower_raw/PlayaSpring_raw_GPF_LinDet_2013_05_24.mat'));
        t_start=3236858+18000*chunk_15;   %written for 15 min chucks at 20 Hz
        t_end = t_start + 18000;
        START_TIME=datetime(datevec(Sonic_data.playaRawData.t(t_start,1))) %double checks starting at right time
        %Define tower vars
        %col 4-2.02m data
        %col 5-0.61m data
        w_prime = Sonic_data.playaRawData.wPF_Prime(t_start:t_end,5); % mean(TowerData.playaRawData.wPF_Prime(tstart:tend,:),2);
        T_prime = Sonic_data.playaRawData.sonTsPrime(t_start:t_end,5);% mean(TowerData.playaRawData.sonTsPrime(tstart:tend,:),2);
        clear rawFlux;
        
        %air properties
        Cp_air = 1.007e3; % [J kg^-1 K^-1]300k Heat Trans book
        rho_air = 1.1614; %[kg m^-3] 300K HT book
end


%Properties and Constants
%camera properties
dt = 1/20; %20 Hz
Nx=size(T_surface,1);%IOP10~256;
Ny=size(T_surface,2);%IOP10~320;
Nt = size(T_surface,3);
xPix = 1:1:Nx;
yPix = 1:1:Ny;

%soil properties
rho_soil = 1426; %? [kg/m^-3] quick look/Eric
Cp_soil = 2.2e6; % [J m^-3 K^-1] 
alpha_soil = .4e-6; % [m^2 s^-1]
k_soil = .9; % [W m^-2 K^-1]
tau_soil=Avg_chunk; %[s] 
dz = sqrt(4*alpha_soil*tau_soil); %[m]

%Fore Restore Variables
T_deepsoil=290.15; %temperature taken at 70cm below surface
Omega = (2*pi)/(24*3600); % angular rotational speed of the earth [rad/sec]
Cg = Cp_soil*sqrt(alpha_soil/(2*Omega));
mu_s= Omega * Cg; 

%% Apply spatial correction for albedo

%calculate TIR components
avg_window = Nt;
[T_fluct, T_patch, T_trend, T_mean] = compute_TIR_components(T_surface,avg_window);
%%
mean_albedo = mean(SW_out(1:3))./mean(SW_in(1:3));
f = 15;
albedo_star = zeros(Nx,Ny);

for x = 1:Nx
    for y = 1:Ny
        if T_patch(x,y) <= 0 
            albedo_star(x,y) = mean_albedo.*(1 - (((T_patch(x,y))/T_mean)*f));
        elseif T_patch(x,y) > 0
            albedo_star(x,y) = mean_albedo.*(1 - (((T_patch(x,y))/T_mean)*f));
        end
    end
end

Net_Radiation = mean(SW_in(1:3)).*(1-albedo_star)+mean(LW_in(1:3))-mean(LW_out(1:3));

figure;pcolor(Net_Radiation);colorbar
figure;pcolor(albedo_star);colorbar

%figure;pcolor(T_patch);colorbar

%% SEB-TIR calculator

[ q_air, q_soil, E_soil ] = SEB_TIR_calc( tau_soil, dt, Cp_soil, dz, mu_s ,T_deepsoil, T_surface, Net_Radiation);


%% Stats
RUNTIME=toc
q_mean=squeeze(mean(mean(q_air,1,'omitnan'),2,'omitnan'));
q_std=squeeze(mean(std(q_air,'omitnan'),'omitnan'));
q_mean_soil=squeeze(mean(mean(q_soil,1,'omitnan'),2,'omitnan'));
q_std_soil=squeeze(mean(std(q_soil,'omitnan'),'omitnan'));
E_std_soil=squeeze(mean(mean(E_soil,1,'omitnan'),2,'omitnan'));
E_mean_soil=squeeze(mean(mean(E_soil,1,'omitnan'),2,'omitnan'));
q_air_prime = abs(q_mean - mean(q_mean));

%% Make various plots and save workspace

% Save Workspace
%save(sprintf('%sPeriod_wTanalysis__AVG%isec',save_path,Avg_chunk),'-v7.3')

% Plots
% wT tower vs SHF from thermal camera
xline = 0:.0001:250;
yline = 0:.0001:250;
figure()
plot(q_mean,rho_air*Cp_air*mean_wTtower,'o')
hold on 
plot(xline,yline)
hold off
xlabel('$q_{air}$ Thermal Camera','interpreter','latex', 'fontsize',18);
ylabel('$w^\prime T^\prime$ Tower','interpreter','latex','fontsize',18);
%axis([1E-1 50 1E-1 50]);
savefig(sprintf('%s_wTtower_vsTC__AVG%isec.fig',save_path,Avg_chunk))

%Time series of camera vs tower
figure()
plot(q_mean)
hold on 
plot(mean_wTtower*rho_air*Cp_air)
hold off
ylabel('$<wT>_{xy}$  $[Wm^-2]$','interpreter','latex','fontsize',18)
xlabel('$time$ $[s]$','interpreter','latex','fontsize',18)
%title('$<wT>_{xy}$','interpreter','latex','fontsize',18')
legend('Thermal Camera','Experimental Data')
savefig(sprintf('%s_wTvsTime_vsTC__AVG%isec.fig',save_path,Avg_chunk))

%Energy Buget
figure()
plot(linspace(0,900,length(q_mean)),q_mean)
hold on 
plot(linspace(0,900,3),Net_Radiation(1:3))
plot(linspace(0,900,length(q_mean)),q_mean_soil)
plot(linspace(0,900,length(q_mean)),E_mean_soil)
hold off
ylabel('$<E>_{xy}$  $[Wm^-2]$','interpreter','latex','fontsize',18)
xlabel('$time$ $[s]$','interpreter','latex','fontsize',18)
%title('$<wT>_{xy}$','interpreter','latex','fontsize',18')
legend('wT','Net Radiation','Ground Heat Flux','Storage')
savefig(sprintf('%s_Energy_Buget_vsTC__AVG%isec.fig',save_path,Avg_chunk))

% wT tower vs SHF from thermal camera with error bars
figure()
errorbar(q_mean(:),q_std(:),'o')
hold on 
plot(mean_wTtower*rho_air*Cp_air)
hold off
ylabel('$<wT>_{xy}$ $[Wm^-2]$','interpreter','latex','fontsize',18)
xlabel('$time$ $[s]$','interpreter','latex','fontsize',18)
legend('Thermal Camera','Experimental Data')
savefig(sprintf('%s_wTxy_vsTC__AVG%isec.fig',save_path,Avg_chunk))


