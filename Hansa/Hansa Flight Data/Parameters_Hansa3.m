%% Aerodynamic derivatives SET 01 : source bala
% LOAD INITIAL CONSTANT PARAMETER OF AIRCRAFT HANSA-3
%%  Aircraft parameter
c      = 1.211;  % Mean aerodynamic chord in mtr
b      = 10.47;  % Wing span in mtr
AR     = 8.8;    % Aspect ratio
mass   = 750;    % in Kg
zcg    = 0.335;  % Distance of CG of Aircraft in verticle plane in mtr
e      = 0.7885; % oswald efficiancy factor
g      = 9.81;   % in m/s^2
rho    = 1.225;  % density

Ixx    = 873;    % In Kg-m-sq
Iyy    = 907;    % In Kg-m-sq
Izz    = 1680;   % In Kg-m-sq
Ixz    = 1144;  % In Kg-m-sq

L0     = Ixx*Izz-Ixz^2;
L1     = ((Ixx-Iyy+Izz)*Ixz)/L0;
L2     = (Izz*(Izz-Iyy)+ Ixz^2)/L0;
L3     = Izz/L0;
L4     = Ixz/L0;
L5     = (Izz-Ixx)/Iyy;
L6     = Ixz/Iyy;
L7     = (Ixx*(Ixx-Iyy)+ Ixz^2)/L0;
L8     = Ixx/L0;

%% Aerodynamic cofficients 
CL0         = 0.354;
CL_alpha    = 4.972;
CL_q        = 37.259;
CL_delta_e  = 0.265;% 0.54
%
CD0         = 0.036;
CD_alpha    = 0.061;
CD_q        = 0 ;
CD_delta_e  = 0.026;
%
Cm0         =  0.052;
Cm_alpha    = -0.496;
Cm_q        = -11.729;
Cm_delta_e  = -1.008;

%
CY0         = -0.013; 
CY_beta     = -0.531;
CY_p        = 0.3269;
CY_r        = 0.633;
CY_delta_a  = 0;
CY_delta_r  = 0.150;
%
Cl0         =  0.0015;
Cl_beta     = -0.031;
Cl_p        = -0.4251;
Cl_r        =  0.1836; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cl_delta_a  = -0.153;
Cl_delta_r  =  0.005;
%
Cn0         =  0.001;
Cn_beta     =  0.061;
Cn_p        = -0.055;
Cn_r        = -0.091;
Cn_delta_a  =  0;
Cn_delta_r  = -0.049;

%% Parameter required for the atmosphere module
T0 = 15;  % Intial Temperature sea level deg celcius
h0 =  0;  % Intial altitude sea level in km
h1 = 11;  % Intial altitude sea level in km
h2 = 20;  % Intial altitude sea level in km
