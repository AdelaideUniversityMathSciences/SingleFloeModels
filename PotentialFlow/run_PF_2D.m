% function [s_s,s_h,s_p] = run_PF_2D(f,h,d,D,L,N,M,rho,A_p0,B_m0,modes,TS)
%
% LJ YIEW
% Created on  Feb 2014
% Last edited Oct 2016
%
% Solves the 2D problem for a floe excited by an incident wave.
% Calculates the surge, heave and pitch amplitudes.
%
% INPUTS:
% f     = frequency [Hz]
% h     = water depth [m]
% d     = draft [m]
% D     = floe thickness [m]
% L     = half length of floe [m]
% N     = number of modes (vertical eigenfunctions)
% M     = mass of floe [kg]
% rho   = density of fluid [kg/m^3]
% A_p0  = incident wave amplitude (REGION 1) [m]
% B_m0  = incident wave amplitude (REGION 3) [m]
% modes = degrees of freedom (surge,heave,pitch) 
%         eg. 100 = surge only, 101 = surge & pitch only
% TS    = flag for higher order thickness terms
%
% OUTPUTS:
% [s_s,s_h,s_p] = complex surge,heave,pitch amplitudes
%
% FILES NEEDED:
%  wavefield.m
%  dispersion.m
%  dispersion_rad.m
%  fn_Diffraction.m
%  fn_Surge.m
%  fn_Heave.m
%  fn_Pitch.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLUID REGIONS:
% 
%          A_p ==>                                 <== B_m
%
%              ----------------------------------------- z=0
%                         |                 |
%                         |      FLOE       |
%                         |_________________|            z=-d
%                         :                 :
%                         :                 :
%                REGION 1 :     REGION 2    : REGION 3
%              ___________:_________________:___________ z=-h
%                  
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s_s,s_h,s_p] = run_PF_2D(f,h,d,D,L,N,M,rho,A_p0,B_m0,modes,TS)

if ~exist('f','var'); f=0.5; end
if ~exist('h','var'); h=0.8; end
if ~exist('d','var'); d=0.01; end
if ~exist('D','var'); D=0.015; end
if ~exist('L','var'); L=0.2; end
if ~exist('N','var'); N=20; end
if ~exist('M','var'); M=1; end
if ~exist('rho','var'); rho=1000; end
if ~exist('A_p0','var'); A_p0=1; end
if ~exist('B_m0','var'); B_m0=0; end
if ~exist('modes','var'); modes=111; end
if ~exist('TS','var'); TS=1; end

% INPUTS
s_u = 1;  % unit heave,surge,pitch amplitude
g = 9.81;
omega = 2*pi*f;
sigma = omega^2/g;
[field] = wavefield('f',f,h);
k = cell2mat(field(5,2));
a_p0 = A_p0*exp(-1i*k*L);
b_m0 = B_m0*exp(-1i*k*L);
a_p = [a_p0;zeros(N,1)];
b_m = [b_m0;zeros(N,1)];

% HYDRODYNAMIC TERMS (INTEGRALS OF DIFFRACTION & RADIATION POTENTIALS)
DF_n1 = int_phi_d_n1(sigma,h,d,L,N,a_p,b_m,A_p0,B_m0);
DF_n3 = int_phi_d_n3(sigma,h,d,L,N,a_p,b_m,A_p0,B_m0);
DF_n5 = int_phi_d_n5(sigma,h,d,D,L,N,a_p,b_m,A_p0,B_m0);
SG_n1 = int_phi_s_n1(sigma,s_u,h,d,L,N);
SG_n3 = int_phi_s_n3(sigma,s_u,h,d,L,N);
SG_n5 = int_phi_s_n5(sigma,s_u,h,d,D,L,N);
HV_n1 = int_phi_h_n1(sigma,s_u,h,d,L,N);
HV_n3 = int_phi_h_n3(sigma,s_u,h,d,L,N);
HV_n5 = int_phi_h_n5(sigma,s_u,h,d,D,L,N);
PT_n1 = int_phi_p_n1(sigma,h,d,D,L,N,s_u);
PT_n3 = int_phi_p_n3(sigma,h,d,D,L,N,s_u);
PT_n5 = int_phi_p_n5(sigma,h,d,D,L,N,s_u);

% SECOND MOMENTS OF INERTIA
I_33 = 1/6*d*D^2*L;
I_33 = I_33*TS;
I_11 = 2/3*L^3*d;

% INTEGRAL OF HYDROSTATIC TERM (PITCH)
U_5  = 2/3*L^3 + TS*L*d*(d-D);

% DEGREES OF FREEDOM
switch modes
 case 111
  % SURGE,HEAVE,PITCH
  LHS = [ -M*sigma/rho-SG_n1  -HV_n1                 -PT_n1 ;
          -SG_n3              2*L-M*sigma/rho-HV_n3  -PT_n3 ;
          -SG_n5              -HV_n5                 -sigma*(I_11+I_33)+U_5-PT_n5]; 
  RHS = [ DF_n1 ;
          DF_n3 ;
          DF_n5 ];
  SOL = LHS\RHS;
  s_s = SOL(1);
  s_h = SOL(2);
  s_p = SOL(3);
  
 case 110 
  % SURGE,HEAVE
  LHS = [ -M*sigma/rho-SG_n1  -HV_n1                 ;
          -SG_n3              2*L-M*sigma/rho-HV_n3  ]; 
  RHS = [ DF_n1 ;
          DF_n3 ];
  SOL = LHS\RHS;
  s_s = SOL(1);
  s_h = SOL(2);
  s_p = 0;
  
 case 101
  % SURGE,PITCH
  LHS = [ -M*sigma/rho-SG_n1  -PT_n1 ;
          -SG_n5              -sigma*(I_11+I_33)+U_5-PT_n5];
  RHS = [ DF_n1  ;
          DF_n5 ];
  SOL = LHS\RHS;
  s_s = SOL(1);
  s_h = 0;
  s_p = SOL(2);
  
 case 011
  % HEAVE,PITCH
  LHS = [ 2*L-M*sigma/rho-HV_n3  -PT_n3 ;
          -HV_n5                 -sigma*(I_11+I_33)+U_5-PT_n5];
  RHS = [ DF_n3  ;
          DF_n5 ];
  SOL = LHS\RHS;
  s_s = 0;
  s_h = SOL(1);
  s_p = SOL(2);
  
 case 100
  % SURGE
  LHS = [ -M*sigma/rho-SG_n1 ];  
  RHS = [ DF_n1 ];
  SOL = LHS\RHS;
  s_s = SOL(1);
  s_h = 0;
  s_p = 0;
  
 case 010
  % HEAVE
  LHS = [ 2*L-M*sigma/rho-HV_n3 ];  
  RHS = [ DF_n3 ];
  SOL = LHS\RHS;
  s_s = 0;
  s_h = SOL(1);
  s_p = 0;
  
 case 001
  % PITCH
  LHS = [ -sigma*(I_11+I_33)+U_5-PT_n5 ];  
  RHS = [ DF_n5 ];
  SOL = LHS\RHS;
  s_s = 0;
  s_h = 0;
  s_p = SOL(1);
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% INTEGRALS OF POTENTIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% =========================================================================
% INTEGRAL OF phi wrt n1 dS
% =========================================================================
% DIFFRACTION
function out = int_phi_d_n1(sigma,h,d,L,N,a_p,b_m,A_p0,B_m0)
clear a_m b_p alpha_p alpha_m A B k kappa int_eta int_phi_dz out
[a_m,b_p,alpha_p,alpha_m,A,B,k,kappa] = fn_Diffraction(sigma,h,d,L,N,A_p0,B_m0);                                        
int_eta = (sinh(k*h)-sinh(k*(h-d)))./(k.*cosh(k*h));
for n = 1:N+1
 int_phi_dz(n) = (a_p(n)+a_m(n)-b_p(n)-b_m(n))*int_eta(n); % CHECK SIGNS
end
out = sum(int_phi_dz);
% =========================================================================
% SURGE
function out = int_phi_s_n1(sigma,s,h,d,L,N)
clear a_m b_p alpha_p alpha_m A B u v k kappa int_eta int_phi_dz out
[a_m,b_p,alpha_p,alpha_m,A,B,k,kappa] = fn_Surge(sigma,h,d,L,N,s);
int_eta = (sinh(k*h)-sinh(k*(h-d)))./(k.*cosh(k*h));
for n = 1:N+1
 int_phi_dz(n) = (a_m(n)-b_p(n))*int_eta(n);
end
out = sum(int_phi_dz);
% =========================================================================
% HEAVE
function out = int_phi_h_n1(sigma,s,h,d,L,N)
clear a_m b_p alpha_p alpha_m A B k kappa lambda int_eta int_phi_dz out
[a_m,b_p,alpha_p,alpha_m,A,B,k,kappa] = fn_Heave(sigma,h,d,L,N,s);
int_eta = (sinh(k*h)-sinh(k*(h-d)))./(k.*cosh(k*h));
for n = 1:N+1
 int_phi_dz(n) = (a_m(n)-b_p(n))*int_eta(n);
end
check = sum(int_phi_dz);
out = 0;% sum(int_phi_dz); % check that out = 0
% =========================================================================
% PITCH
function out = int_phi_p_n1(sigma,h,d,D,L,N,s_p)
clear a_m b_p alpha_p alpha_m A B k kappa lambda int_eta int_phi_dz out
[a_m,b_p,alpha_p,alpha_m,A,B,k,kappa] = fn_Pitch(sigma,h,d,D,L,N,s_p);
int_eta = (sinh(k*h)-sinh(k*(h-d)))./(k.*cosh(k*h));
for n = 1:N+1
 int_phi_dz(n) = (a_m(n)-b_p(n))*int_eta(n);
end
out = sum(int_phi_dz);
% =========================================================================


% =========================================================================
% INTEGRAL OF phi wrt n3 dS
% =========================================================================
% DIFFRACTION
function out = int_phi_d_n3(sigma,h,d,L,N,a_p,b_m,A_p0,B_m0)
clear a_m b_p alpha_p alpha_m A B k kappa k_hat kappa_hat int_phi_dx out
[a_m,b_p,alpha_p,alpha_m,A,B,k,kappa] = fn_Diffraction(sigma,h,d,L,N,A_p0,B_m0);

for n = 2:N+1
 int_phi_dx(n) = (alpha_p(n)+alpha_m(n))*(exp(2i*kappa(n)*L)-1)/(1i*kappa(n));
end
out = 2*B*L + sum(int_phi_dx);
% =========================================================================
% SURGE
function out = int_phi_s_n3(sigma,s,h,d,L,N)
clear a_m b_p alpha_p alpha_m A B u v k kappa int_phi_dx out
[a_m,b_p,alpha_p,alpha_m,A,B,k,kappa] = fn_Surge(sigma,h,d,L,N,s);
for n = 2:N+1
 int_phi_dx(n) = (alpha_p(n)+alpha_m(n))*(exp(2i*kappa(n)*L)-1)/(1i*kappa(n));
end
out = 2*B*L;% +  sum(int_phi_dx); % check out = 2*B*L
% =========================================================================
% HEAVE
function out = int_phi_h_n3(sigma,s,h,d,L,N)
clear a_m b_p alpha_p alpha_m A B k kappa lambda int_phi_dx out
[a_m,b_p,alpha_p,alpha_m,A,B,k,kappa] = fn_Heave(sigma,h,d,L,N,s);
for n = 2:N+1
 int_phi_dx(n) = (alpha_p(n)+alpha_m(n))*(exp(2i*kappa(n)*L)-1)/(1i*kappa(n));
end
out = 2*B*L +  sum(int_phi_dx) + sigma/(2*(h-d))*(2*L*(h-d)^2-2/3*L^3);
% =========================================================================
% PITCH
function out = int_phi_p_n3(sigma,h,d,D,L,N,s_p)
clear a_m b_p alpha_p alpha_m A B k kappa lambda int_phi_dx out
[a_m,b_p,alpha_p,alpha_m,A,B,k,kappa] = fn_Pitch(sigma,h,d,D,L,N,s_p);
for n = 2:N+1
 int_phi_dx(n) = (alpha_p(n)+alpha_m(n))*(exp(2i*kappa(n)*L)-1)/(1i*kappa(n));
end
% out = 2*B*L +  sum(int_phi_dx);
out = 2*B*L;
% =========================================================================


% =========================================================================
% INTEGRAL OF phi wrt n5 dS
% =========================================================================
% DIFFRACTION
function out = int_phi_d_n5(sigma,h,d,D,L,N,a_p,b_m,A_p0,B_m0)
clear a_m b_p alpha_p alpha_m A B k kappa k_hat kappa_hat ...
 int_phi_dx int_phi_dz int_eta out
[a_m,b_p,alpha_p,alpha_m,A,B,k,kappa] = fn_Diffraction(sigma,h,d,L,N,A_p0,B_m0);

for n = 2:N+1
 int_phi_dx(n) = (alpha_p(n)-alpha_m(n))...
                   *(L/(1i*kappa(n))*(exp(2i*kappa(n)*L)+1)...
                       +1/kappa(n)^2*(exp(2i*kappa(n)*L)-1));
end
int_eta = ((2*d-D).*sinh(k.*h)+D.*sinh(k.*(h-d)))./(2.*k.*cosh(k.*h))...
             - (cosh(k.*h)-cosh(k.*(h-d)))./(k.^2.*cosh(k.*h));
for n = 1:N+1
 int_phi_dz(n) = (a_p(n)+a_m(n)-b_p(n)-b_m(n))*int_eta(n); % #####
end
out = -(2/3*A*L^3 + sum(int_phi_dx)) + sum(int_phi_dz);
% test = 2/3*A*L^3 + sum(int_phi_dx)
% =========================================================================
% SURGE
function out = int_phi_s_n5(sigma,s,h,d,D,L,N)
clear a_m b_p alpha_p alpha_m A B u v k kappa ...
 int_eta int_phi_dx int_phi_dz out
[a_m,b_p,alpha_p,alpha_m,A,B,k,kappa] = fn_Surge(sigma,h,d,L,N,s);
for n = 2:N+1
 int_phi_dx(n) = (alpha_p(n)-alpha_m(n))...
                   *(L/(1i*kappa(n))*(exp(2i*kappa(n)*L)+1)...
                       +1/kappa(n)^2*(exp(2i*kappa(n)*L)-1));
end
int_eta = ((2*d-D).*sinh(k.*h)+D.*sinh(k.*(h-d)))./(2.*k.*cosh(k.*h))...
             - (cosh(k.*h)-cosh(k.*(h-d)))./(k.^2.*cosh(k.*h));
for n = 1:N+1
 int_phi_dz(n) = (a_m(n)-b_p(n))*int_eta(n);
end
out = -(2/3*A*L^3 + sum(int_phi_dx)) + sum(int_phi_dz);
% =========================================================================
% HEAVE
function out = int_phi_h_n5(sigma,s,h,d,D,L,N)
clear a_m b_p alpha_p alpha_m A B k kappa lambda ...
 int_eta int_phi_dx int_phi_dz out
[a_m,b_p,alpha_p,alpha_m,A,B,k,kappa] = fn_Heave(sigma,h,d,L,N,s);
for n = 2:N+1
 int_phi_dx(n) = (alpha_p(n)-alpha_m(n))...
                   *(L/(1i*kappa(n))*(exp(2i*kappa(n)*L)+1)...
                       +1/kappa(n)^2*(exp(2i*kappa(n)*L)-1));
end
int_eta = ((2*d-D).*sinh(k.*h)+D.*sinh(k.*(h-d)))./(2.*k.*cosh(k.*h))...
             - (cosh(k.*h)-cosh(k.*(h-d)))./(k.^2.*cosh(k.*h));
for n = 1:N+1
 int_phi_dz(n) = (a_m(n)-b_p(n))*int_eta(n);
end
% check = sum(int_phi_dz)
% out = 2/3*A*L^3 + sum(int_phi_dx)% + sum(int_phi_dz); % check sum(int_phi_dz) = 0
out = -(2/3*A*L^3);
% =========================================================================
% PITCH
function out = int_phi_p_n5(sigma,h,d,D,L,N,s_p)
clear a_m b_p alpha_p alpha_m A B k kappa lambda ...
 int_eta int_phi_dx int_phi_dz out
[a_m,b_p,alpha_p,alpha_m,A,B,k,kappa] = fn_Pitch(sigma,h,d,D,L,N,s_p);
for n = 2:N+1
 int_phi_dx(n) = (alpha_p(n)-alpha_m(n))...
                   *(L/(1i*kappa(n))*(exp(2i*kappa(n)*L)+1)...
                       +1/kappa(n)^2*(exp(2i*kappa(n)*L)-1));
end
int_eta = ((2*d-D).*sinh(k.*h)+D.*sinh(k.*(h-d)))./(2.*k.*cosh(k.*h))...
             - (cosh(k.*h)-cosh(k.*(h-d)))./(k.^2.*cosh(k.*h));
for n = 1:N+1
 int_phi_dz(n) = (a_m(n)-b_p(n))*int_eta(n); % ###
end
out = -(2/3*A*L^3 + sum(int_phi_dx) - sigma*1*((h-d)*L^3/3 - L^5/(15*(h-d))))...
      + sum(int_phi_dz);
% =========================================================================