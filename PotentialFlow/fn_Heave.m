% function [a_m,b_p,alpha_p,alpha_m,A,B,k,kappa] = fn_Heave(sigma,h,d,L,N,s)
%
% LJ YIEW
% Created on  Nov 2013
% Last edited Oct 2016
%
% Solves the radiation problem for a floe excited by heave, calculates the
% radiated wave coefficients.
%
% INPUTS:
% sigma = frequency parameter omega^2/g where omega = 2*pi*f
% h     = water depth [m]
% d     = draft [m]
% L     = half length of floe [m]
% N     = number of modes (vertical eigenfunctions)
% s     = heave amplitude [m]
%
% OUTPUTS:
% a_m   = reflected wave coefficients (REGION 1)
% b_p   = transmitted wave coefficients (REGION 1)
% A,B,alpha_p,alpha_m = wave coefficients (REGION 2)
% k     = wave numbers (REGION 1&3)
% kappa = wave numbers (REGION 2)
%
% NB.
% REGION 1&3 ARE THE FREE SURFACE REGIONS (TO THE LEFT & RIGHT OF THE FLOE)
% REGION 2 IS THE REGION UNDER THE FLOE
%
% FILES NEEDED:
%  dispersion_rad.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a_m,b_p,alpha_p,alpha_m,A,B,k,kappa] = fn_Heave(sigma,h,d,L,N,s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLUTION MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ROOTS OF DISPERSION EQUATION
[k,kappa] = dispersion_rad(sigma,h,d,N,L);

% MATRIX OF INNER PRODUCTS
for m = 1:N+1
 Om(m,1) = O_m(kappa(m),h,d);
 F1(m,1) = F_m(h,d,L,kappa(m))*sigma*s/(2*(h-d));
 F2(m,1) = O_m(kappa(m),h,d)*sigma*s*L/(h-d);
 for n = 1:N+1
  O(m,n) = O_mn(kappa(n),kappa(m),h,d);
  D(m,n) = D_mn(k(n),kappa(m),h,d);
  C(m,n) = D_mn(k(m),kappa(n),h,d);
  Q(m,n) = Q_mn(k(m),k(n),h);
 end
end
O2 = O(:,2:end);
% OTHER FUNCTIONS
K = diag(1i*k);
P = diag(1i*kappa(2:end));
E = diag(exp(2i*kappa(2:end)*L));
% ZERO MATRICES
zero  = zeros(N+1);
zero1 = zeros(N+1,1);
zeroN = zeros(N+1,N);

% COEFFICIENTS OF RADIATION TERMS
S = [ D     zero  Om*L   -O2     -Om    -O2*E    zero  zero ;
      -Q*K  zero  zero1  zeroN   zero1  zeroN    -C    zero ;
      zero  zero  Om     O2*P    zero1  -O2*P*E  -O    zero ;
      zero  D     -Om*L  -O2*E   -Om    -O2      zero  zero ;
      zero  Q*K   zero1  zeroN   zero1  zeroN    zero  -C   ;
      zero  zero  Om     O2*P*E  zero1  -O2*P    zero  -O   ];

% FORCING MATRIX
 I = [ F1    ;
       zero1 ;
       -F2   ;
       F1    ;
       zero1 ;
       F2    ];

% SOLVE
R = inv(S)*I;

% EXTRACT RADIATION SOLUTION
a_m     = R(1:N+1);
b_p     = R(N+1+1:2*(N+1));
A       = R(2*(N+1)+1);
alpha_p = R(2*(N+1)+2:3*(N+1));
B       = R(3*(N+1)+1);
alpha_m = R(3*(N+1)+2:4*(N+1));
u       = R(4*(N+1)+1:5*(N+1));
v       = R(5*(N+1)+1:6*(N+1));
alpha_p = [0;alpha_p];
alpha_m = [0;alpha_m];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INNER PRODUCTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = Q_mn(k_m,k_n,h)
% syms k h z; int(cosh(k*(z+h))^2/cosh(k*h)^2,z,-h,0)
if k_n == k_m
 out = 1/2*(sinh(k_n*h)*cosh(k_n*h)+k_n*h)/(k_n*cosh(k_n*h)^2);
else
 out = 0;
end

function out = D_mn(k,kappa,h,d)
% syms k kappa h d z; 
% int(cosh(k*(z+h))/cosh(k*h)*cosh(kappa*(z+h))/cosh(kappa*(h-d)),z,-h,-d)
out = -(k*cosh(kappa*(d - h))*sinh(k*(d - h))...
      - kappa*cosh(k*(d - h))*sinh(kappa*(d - h)))...
       /(cosh(h*k)*cosh(kappa*(d - h))*(k^2 - kappa^2));
      
function out = O_mn(kappa_n,kappa_m,h,d)
% syms kappa h d z; int(cosh(kappa*(z+h))^2/cosh(kappa*(h-d))^2,z,-h,-d)
if kappa_n == 0 && kappa_m == 0
 out = h-d;
elseif kappa_n == kappa_m
 out = -(d/2 - h/2 + sinh(2*kappa_n*(d - h))/(4*kappa_n))/cosh(kappa_n*(d - h))^2;
else
 out = 0;
end

function out = O_m(kappa,h,d)
% syms kappa z h d; out = int(cosh(kappa*(z+h))/cosh(kappa*(h-d)),z,-h,-d)
if kappa == 0
 out = h-d;
else
 out = 1/kappa*tanh(kappa*(h-d));
end

function out = F_m(h,d,L,kappa)
% syms z h L kappa d;
% int(((z+h)^2-L^2)*cosh(kappa*(z+h))/cosh(kappa*(h-d)),z,-h,-d)
% int((z+h)^2-L^2,z,-h,-d)
if kappa == 0
 out = ((d - h)*(3*L^2 - d^2 + 2*d*h - h^2))/3;
else
 out = (sinh(kappa*(d - h))*L^2)/(kappa*cosh(kappa*(d - h)))...
       - (2*sinh(kappa*(d - h)) - 2*kappa*cosh(kappa*(d - h))*(d - h)...
          + kappa^2*sinh(kappa*(d - h))*(d - h)^2)...
          /(kappa^3*cosh(kappa*(d - h)));
end
 


 