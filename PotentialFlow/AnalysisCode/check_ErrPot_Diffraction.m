% function check_ErrPot_Diffraction
%
% LJ YIEW
% Created on  Nov 2014
% Last edited Dec 2016
%
% Considers the diffraction problem.
% Checks the accuracy of the eigenfunction matching method by calculating
% the error between velocities at the interfacial fluid boundaries over a 
% range of wave frequencies and number of modes.
%
% FILES NEEDED:
%  Param_AMC.m
%  fn_Diffraction.m
%  ColourCode.m
%  dispersion.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function check_ErrPot_Diffraction

tic
close all
clear all
clc

addpath('../')

% INPUTS:
Param = Param_AMC;
L     = Param.L; % half length of floe [m]
h     = Param.h; % water depth [m]
d     = Param.d; % draft [m]
D     = Param.D; % thickness [m]

% RANGE OF FREQUENCY
Freq  = linspace(0.5,2,50);
Omega = Freq*2*pi;
g     = 9.81;
Sigma = Omega.^2/g;
% RANGE OF Z
Z = linspace(-h,-d,2000);
Z_interval = diff(Z);
% RANGE OF N
N_range = 20:5:400;

%%
% LOOP THRU FREQUENCY
for j = 1:length(Omega)
 
 freq = Freq(j);
 sigma = Sigma(j);
 
 % LOOP THRU N
 for jj = 1:length(N_range)
  
  N = N_range(jj)
  
  % INCIDENT WAVE AMPLITUDE
  A_p0  = 1;   % amplitude of incident wave in +ve x direction
  A_p   = [A_p0;zeros(N,1)];
  B_m0  = 0;   % amplitude of incident wave in -ve x direction
  B_m   = [B_m0;zeros(N,1)];
  % add phase shift
  [k,~] = dispersion(sigma,h,d,N);
  a_p   = diag(exp(-1i*k*L))*A_p;
  b_m   = diag(exp(-1i*k*L))*B_m;

  % OUTPUTS:
  clear a_m b_p alpha_p alpha_m A B k kappa
  [a_m,b_p,alpha_p,alpha_m,A,B,k,kappa] = fn_Diffraction(sigma,h,d,L,N,A_p0,B_m0);
  
  % LOOP THRU Z
  for jjj = 1:length(Z)
   
   z = Z(jjj);
   
   % POTENTIAL IN REGION 1
   for n=1:N+1
    phi(n) = (a_p(n)+a_m(n))*cosh(k(n)*(z+h))/cosh(k(n)*h);
   end
   phi_1 = abs(sum(phi));
   clear phi
   % POTENTIAL IN REGION 2 (x = -L)
   for n = 2:N+1
    phi(n) = (alpha_p(n)+alpha_m(n)*exp(2i*kappa(n)*L))...
              *cosh(kappa(n)*(z+h))/cosh(kappa(n)*(h-d));
   end
   phi_2L = abs(sum(phi)+A*(-L)+B);
   clear phi
   % POTENTIAL IN REGION 2 (x = +L)
   for n = 2:N+1
    phi(n) = (alpha_p(n)*exp(2i*kappa(n)*L)+alpha_m(n))...
              *cosh(kappa(n)*(z+h))/cosh(kappa(n)*(h-d));
   end
   phi_2R = abs(sum(phi)+A*L+B);
   clear phi
   % POTENTIAL IN REGION 3
   for n=1:N+1
    phi(n) = (b_p(n)+b_m(n))*cosh(k(n)*(z+h))/cosh(k(n)*h);
   end
   phi_3 = abs(sum(phi));
   clear phi
   
   % ERROR
   errorL(jjj) = abs(phi_1 - phi_2L);
   errorR(jjj) = abs(phi_3 - phi_2R);
   clear phi_1 phi_2L phi_2R phi_3
   
  end
  
  EL(jj,j) = sum(errorL)*Z_interval(1);
  ER(jj,j) = sum(errorR)*Z_interval(1);
  clear errorL errorR
  
 end
 
end

toc

%% 

% PLOT ERROR BETWEEN POTENTIALS VS N

% COLOR RANGE
[colour] = ColourCode(Freq);

for j = 1:length(Freq) 
 
 freq = Freq(j);
 
 figure(1)
 hold on
 plot(0,0,'g',0,0,'b') % set the colours for the legend
 plot(N_range,EL(:,j),'Color',cell2mat(colour(j)))
 set(gcf,'position',[100 100 500 300]);
 set(gca,'FontSize',14)
 set(gca,'yscale','log')
 xlabel('N'),ylabel('log_{10}(Error)')
 legend([num2str(Freq(1)),' Hz'],[num2str(Freq(end)),' Hz'])
 title(['x = -L '])

 figure(2)
 hold on
 plot(0,0,'g',0,0,'b') % set the colours for the legend
 plot(N_range,ER(:,j),'Color',cell2mat(colour(j)))
 set(gcf,'position',[100 100 500 300]);
 set(gca,'FontSize',14)
 set(gca,'yscale','log')
 xlabel('N'),ylabel('log_{10}(Error)')
 legend([num2str(Freq(1)),' Hz'],[num2str(Freq(end)),' Hz'])
 title(['x = +L '])
 
end