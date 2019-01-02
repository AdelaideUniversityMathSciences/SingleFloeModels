% function plot_RAO
%
% LJ YIEW
% Created on  Aug 2017
%
% Generates RAOs for AMC experiments
%
% INPUTS:
%
% OUTPUTS:
% surge,heave,pitch RAOs
%
% FILES NEEDED:
%  wavefield.m
%  run_PF_2D.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_RAO

% % PHYSICAL PARAMETERS FOR AMC EXPERIMENTS
% Lambda = linspace(0.4,6,50);
% h = 0.831;
% L = 0.2;
% rho_b = 650;
% rho = 1000;
% D = 0.015;

% PHYSICAL PARAMETERS FOR BGO (OCEANIDE) EXPERIMENTS
dat = open('param_bgo.mat');
Lambda = dat.bgo_wavelengths;
h = 3.1;
L = 0.495;
rho_b = 551;
rho = 1000;
D = 0.033;

N = 100;

% ARCHIMEDES CONDITION
d = rho_b/rho*D;
m_3D = rho_b*D*pi*0.2^2; % [kg]
m = rho_b*2*L*D;         % [kg/m]

% WAVE PARAMETERS
A_p0  = 1;
B_m0  = 0;

% MOTIONS
modes = 111; % [1,0] CORRESPONDING TO SURGE,HEAVE,PITCH

TS = 1; % THICKNESS TERMS

for j = 1:length(Lambda)
  lambda  = Lambda(j)
  [field] = wavefield('lambda',lambda,h);
  f       = cell2mat(field(1,2));
  k(j)    = cell2mat(field(5,2));
  sigma   = (2*pi*f)^2/9.81;

  [s_s(j),s_h(j),s_p(j)] = run_PF_2D(f,h,d,D,L,N,m,rho,A_p0,B_m0,modes,TS);
 
end

rao_s_2dr = abs(s_s)./A_p0.*tanh(k*h);
rao_h_2dr = abs(s_h)/A_p0;
rao_p_2dr = abs(s_p)./k/A_p0;

%%


figure(1)
hold on
set(gcf,'position',[100 400 500 400]);
set(gca,'FontSize',16)
plot(Lambda/2/L,rao_s_2dr,'ko')
ylabel('Surge RAO')
xlabel('Wavelength / Floe length [m/m]')
ylim([0 1.2])
box on
figure(2)
hold on
set(gcf,'position',[600 400 500 400]);
set(gca,'FontSize',16)
plot(Lambda/2/L,rao_h_2dr,'ko') 
ylabel('Heave RAO')
xlabel('Wavelength / Floe length [m/m]')
ylim([0 1.2])
box on
figure(3)
hold on
set(gcf,'position',[350 0 500 400]);
set(gca,'FontSize',16)
plot(Lambda/2/L,rao_p_2dr,'ko')
ylabel('Pitch RAO')
xlabel('Wavelength / Floe length [m/m]')
ylim([0 1.2])
box on

pause
