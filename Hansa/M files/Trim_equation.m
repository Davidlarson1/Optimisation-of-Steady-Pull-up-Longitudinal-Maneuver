 function Xdot = Trim_equation(X0,T,de)
 
 
 V  = X0(1);
 alpha = X0(2); 
 q = X0(3);
 theta = X0(4);
 x = X0(5);
 z = X0(6);

% Aerodynamic and Geometric Parameters
  GeoPara = [750,9.81 ,1.225,12.47, 1.21, 0.16];
  mass = GeoPara(1);g = GeoPara(2);rho = GeoPara(3);S = GeoPara(4);c_bar = GeoPara(5);
  Iyy = 907; % Remenber to change this while using F-18
 
  AeroPara = [0.354,0.036,0.052,4.972 * 0.0175 , 37.259 * 0.0175 ,-1.008 * 0.0175,- 0.496 * 0.0175,-11.729 * 0.0175,0.265 * 0.0175,0.026 * 0.0175,0.061* 0.0175 ];

CL0 = AeroPara(1);CD0 = AeroPara(2);Cm0 = AeroPara(3);
CL_alpha = AeroPara(4); CL_q = AeroPara(5);
Cm_de = AeroPara(6); Cm_alpha = AeroPara(7);Cm_q = AeroPara(8);CL_de = AeroPara(9);CD_de = AeroPara(10);CD_alpha = AeroPara(11);

  
  

CL = (CL0 + (CL_alpha * alpha) + (CL_q *((q * c_bar)/(2*V))) + CL_de * de);
CD = (CD0 + (CD_alpha * alpha) + CD_de * de);
Cm = Cm0 + (Cm_alpha * alpha) + (Cm_q *((q * c_bar)/(2*V))) + Cm_de * de;



L = 0.5 * rho * V^2 * S * CL;
D = 0.5 * rho * V^2 * S * CD;
M = 0.5 * rho * V^2 * S * Cm * c_bar;

V_dot = -D/mass - g * sind(theta - alpha) + (T * cosd(alpha))/mass;
alpha_dot = -L/(mass * V) + (g * cosd(theta - alpha))/V + q - (T * sind(alpha))/(mass * V);
q_dot  = M/ Iyy;
theta_dot = q;

xdot = (V * cosd(alpha) * cosd(theta)) + (V * sind(alpha) * sind(theta));
zdot = -(V * cosd(alpha) * sind(theta)) + (V * sind(alpha) * cosd(theta));

Xdot = [V_dot,alpha_dot,q_dot,theta_dot,xdot,zdot]';
 
 end