% function run_SS_1Floe(H,lambda)
%
% LJ YIEW & MH MEYLAN
% Created on  Jul 2013
% Last edited Oct 2016
%
% Solves the Rumer/Marchenko slope-sliding model.
%
% INPUTS:
% H      = wave height [m]
% lambda = wavelength [m]
%
% OUTPUTS:
% Plots of x-displacement vs time
%
% FILES NEEDED:
% ParamDef.m
% SS_ode.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_SS_1Floe(H,lambda)

 close all
 clear all
 clc
 
 tic
 
 warning('off')

 if ~exist('H','var'); H = 0.04; end
 if ~exist('lambda','var'); lambda=5*0.4; end

 % DEFINE WAVE/FLOE PARAMETERS
 Param = ParamDef;
 r     = Param.L;          % radius
 h     = Param.h;          % water depth
 dr    = Param.d;          % draft
 D     = Param.D;          % thickness
 g     = Param.g;          % gravity
 A     = pi*r^2+2*pi*r*dr; % wetted surface area
 rho_b = 650;              % floe density
 rho   = 1000;             % fluid density
 m     = rho_b*D*pi*r^2;   % floe mass [kg]

 % SET WAVE PARAMETERS
 k               = 2*pi/lambda;
 omega           = sqrt(g*k*tanh(k*h));
 WaveParam.H     = H;      % wave height
 a               = H/2;    % wave amplitude
 WaveParam.omega = omega;  % angular frequency
 WaveParam.k     = k;      % wavenumber
 WaveParam.rho   = rho;    % fluid density
 WaveParam.h     = h;      % water depth
 
 % SET FLOE PARAMETERS
 FloeParam.m = m;
 FloeParam.A = A;
 
 % SET DAMPING & ADDED MASS
 Cm       = 0.1;
 Cd       = 0;
 Coeff.Cd = Cd;
 Coeff.Cm = Cm; 
 
 % SET MOORING PARAMETERS
 Mooring.K = 0;
 Mooring.C = 0;
 
 % SET TRANSIENT WAVE PARAMETERS
%  f.a = 0;
%  f.b = 0;
%  f.c = 0;
 Trans.f = 0;
 Trans.t = 0; % flag for transient solution (time dependent amplitude)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SINGLE FLOE MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 % SIMULATION PARAMETERS
 tspan = linspace(0,60,1000); % RANGE OF TIME
 IC1   = [0,0]; % INITIAL CONDITIONS [DISPLACEMENT, VELOCITY]

 %% SOLVE NUMERICAL SOLUTION

 [t1,X1] = ode45(@(t,X) ...
             fn_SS_ode(t,X,0,WaveParam,FloeParam,Coeff,Mooring,Trans), ...
             tspan,IC1);
              
 %% 
 % PLOT X VS TIME
 figure
 set(gcf,'position',[100 400 1000 400]);
 set(gca,'FontSize',14)
 hold on
 plot(t1,X1(:,1)*1e3,'b-') % DISPLACEMENT
 title(['Wavelength = ',num2str(lambda),' m, Wave Amplitude = ',num2str(a),' m'])
 ylabel('x [mm]')
 xlabel('t [s]')
 ylim('auto')
 box on
 grid on
 hold off

toc

end
