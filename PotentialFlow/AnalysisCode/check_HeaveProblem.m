clear all
clc

addpath('../')

f = 0.3;
h = 5;
d = 1;
D = 1;
L = 1;
N = 20;
rho = 1000;
A_p0 = 1;
B_m0 = 0;
modes = 010;
TS = 1;

M(1) = 1;
M(2) = 100;


for j = 1:2
 MM = M(j); 
[s_s,s_h,s_p] = run_PF_2D(f,h,d,D,L,N,MM,rho,A_p0,B_m0,modes,TS);

sh(j) = s_h;
  
 sigma = (2*pi*f)^2/9.81;

  % REFLECTION & TRANSMISSION COEFFICIENTS:
  [a_m_d,b_p_d,~,~,~,~,~,~] = fn_Diffraction(sigma,h,d,L,N,A_p0,B_m0);                                                     
  [a_m_s,b_p_s,~,~,~,~,~,~] = fn_Surge(sigma,h,d,L,N,s_s);
  [a_m_h,b_p_h,~,~,~,~,~,~] = fn_Heave(sigma,h,d,L,N,s_h);
  [a_m_p,b_p_p,~,~,~,~,~,~] = fn_Pitch(sigma,h,d,D,L,N,s_p);
  TC(j) = (b_p_d(1) + b_p_s(1) + b_p_h(1) + b_p_p(1));
  RC(j) = (a_m_d(1) + a_m_s(1) + a_m_h(1) + a_m_p(1));
  
  ad(j) = a_m_d(1);
  
end
  
 %%
 
 LHS = (RC(1) - ad(1))/(RC(2) - ad(2))
 
 RHS = sh(1)/sh(2)
 
 
 
 
 
 
  
  
  
