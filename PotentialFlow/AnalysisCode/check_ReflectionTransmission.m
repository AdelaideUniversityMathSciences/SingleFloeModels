

function check_ReflectionTransmission(lambda)

addpath('../')

% SIMULATION PARAMETERS:
Param = Param_AMC;
L     = Param.L; % half length of floe [m]
h     = Param.h; % water depth [m]
d     = Param.d; % draft [m]
D     = Param.D; % thickness [m]
rho   = Param.rho; % fluid density [kg/m^3]
M     = Param.M; % [kg/m]
N     = 20;

% WAVE PARAMETERS:
if ~exist('lambda','var'); lambda = 1; end
[field] = wavefield('lambda',lambda,h);
f       = cell2mat(field(1,2));
k       = cell2mat(field(5,2));
sigma   = (2*pi*f)^2/9.81;

% INCIDENT WAVE AMPLITUDES:
A_p0 = 1;
B_m0 = 0;


% 2D PF MODEL:
[s_s,s_h,s_p] = run_PF_2D(f,h,d,D,L,N,M,rho,A_p0,B_m0,111,1);
% REFLECTION & TRANSMISSION COEFFICIENTS:
[a_m_d,b_p_d,~,~,~,~,~,~] = fn_Diffraction(sigma,h,d,L,N,A_p0,B_m0);                                                     
[a_m_s,b_p_s,~,~,~,~,~,~] = fn_Surge(sigma,h,d,L,N,s_s);
[a_m_h,b_p_h,~,~,~,~,~,~] = fn_Heave(sigma,h,d,L,N,s_h);
[a_m_p,b_p_p,~,~,~,~,~,~] = fn_Pitch(sigma,h,d,D,L,N,s_p);
TC = abs(b_p_d(1) + b_p_s(1) + b_p_h(1) + b_p_p(1))
RC = abs(a_m_d(1) + a_m_s(1) + a_m_h(1) + a_m_p(1))

checkSolution = TC^2 + RC^2

return