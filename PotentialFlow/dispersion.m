% function [k,kappa] = dispersion(sigma,h,d,N)
%
% LJ YIEW
% Created on  Nov 2013
% Last edited Oct 2016
%
% Finds the real and imaginary roots of open water and ice covered regions
% for the diffraction and surge problem.
%
% INPUTS:
% sigma = frequency parameter omega^2/g where omega = 2*pi*f
% h     = water depth [m]
% d     = draft [m]
% N     = number of modes (vertical eigenfunctions)
%
% OUTPUTS:
% k     = real and imaginary roots of dispersion equation in open water
%         regions
% kappa = imaginary roots in ice covered regions
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [k,kappa] = dispersion(sigma,h,d,N)

% for free surface regions
for n = 1:N
 k(n+1,1) = fzero(@(k) k*tan(k*h)+sigma,...
                       [(2*n-1)/2*pi/h+1e-6,n*pi/h]); % in this interval
 k(n+1,1) = k(n+1,1)*1i;
end
k(1,1) = abs(fzero(@(k) k*tanh(k*h)-sigma,0));

% for dock-covered regions 
% (no-flow condition at z=-d)
for n = 0:N
 kappa(n+1,1) = n*pi/(h-d)*1i;
end

