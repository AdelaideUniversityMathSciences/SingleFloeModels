% function PlotDrift
%
% LJ YIEW
% Created on  Oct 2015
% Last edited Aug 2016
%
% Plots the drift velocities from AMC single-floe experiments
%
% FILES NEEDED:
%  AMC_DataDrift.xlsx
%  wavefield.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotDrift

clear all
clc
% close all


L = 0.2;
h = 0.831;

% READ EXPERIMENTAL DATA
dat_ice = xlsread('AMC_DataDrift.xlsx','A5:L67'); % conditions 1,2,3,4,5 only

% ARRANGE DATA
celerity    = dat_ice(:,7)*1e3; 
driftB      = dat_ice(:,9);
driftNB     = dat_ice(:,11);
steepness   = dat_ice(:,3);
wavelength  = dat_ice(:,6);       
wavenum     = 2*pi./wavelength;
waveheight  = dat_ice(:,5)*1e-3;           % meters
waveheightR = round(waveheight.*1e2)./1e2; % round waveheight to nearest 10mm
% SORT WAVE HEIGHTS
for j = 1:length(waveheightR)
 if waveheightR(j) <= 0.025
  whR(j) = 1;
 elseif waveheightR(j) > 0.025 && waveheightR(j) <= 0.04
  whR(j) = 2;
 elseif waveheightR(j) > 0.04
  whR(j) = 3;
 end
end

%% 
% PLOT DRIFT VELOCITY / PHASE VELOCITY VS WAVE STEEPNESS
% hold on
figure(2)
hold on

plot(steepness,driftB./celerity,'kx','MarkerSize',10)
plot(steepness,driftNB./celerity,'ro','MarkerSize',10)

xlabel('ka')
ylabel('Drift Velocity / Celerity')
box on

% PLOT STOKES DRIFT
k = linspace(0,0.35,10000);
for j = 1:length(k)
 [field]  = wavefield('k',k(j),0.831);
 omega(j) = cell2mat(field(2,2));
 freq(j) = field{1,2};
 lambda(j) = field{4,2};
end

hold on
plot(k,omega.*k./freq./lambda,'k--')







