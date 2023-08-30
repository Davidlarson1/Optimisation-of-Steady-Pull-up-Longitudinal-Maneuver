
 function Xdot = Trim_equations_F18(X0,T,de,AeroPara0)
 
 
 V  = X0(1);
 alpha = X0(2); 
 q = X0(3);
 theta = X0(4);
 x = X0(5);
 z = X0(6);

% Aerodynamic and Geometric Parameters
 GeoPara = [1128.09*14.593 , 9.81 , 1.225 , 400*(0.305^2) , 11.52*0.305 , 239720.815 , 239720.815 ,259969.9570048 ,-2.6, 37.42*0.305];
   mass = GeoPara(1);g = GeoPara(2);rho = GeoPara(3);S = GeoPara(4);c_bar = GeoPara(5);

  Iyy = 239720; % Remenber to change this while using F-18
 
  

CL0 = AeroPara0(5);CD0= AeroPara0(9);Cm0 = AeroPara0(1);
Cm_alpha = AeroPara0(2);CL_alpha = AeroPara0(6);CD_alpha = AeroPara0(10);
CL_q= AeroPara0(7) ;CD_q = AeroPara0(11);Cm_q = AeroPara0(3);
Cm_del = AeroPara0(4); CL_del = AeroPara0(8);CD_del = AeroPara0(12);

AeroPara0 = [Cm0,Cm_alpha,Cm_q,Cm_del,CL0,CL_alpha ,CL_q,CL_del,CD0,CD_alpha,CD_q,CD_del];

  
  

CL = (CL0 + (CL_alpha * alpha) + (CL_q *((q * c_bar)/(2*V))) + CL_del * de);
CD = (CD0 + (CD_alpha * alpha) + CD_del * de + (CD_q *((q * c_bar)/(2*V)))); 
Cm = Cm0 + (Cm_alpha * alpha) + (Cm_q *((q * c_bar)/(2*V))) + Cm_del * de;



L = 0.5 * rho * V^2 * S * CL;
D = 0.5 * rho * V^2 * S * CD;
M = 0.5 * rho * V^2 * S * Cm * c_bar;

V_dot = -D/mass - g * sin(theta - alpha) + (T * cos(alpha))/mass;
alpha_dot = -L/(mass * V) + (g * cos(theta - alpha))/V + q - (T * sin(alpha))/(mass * V);
q_dot  = M/ Iyy;
theta_dot = q;

xdot = (V * cos(alpha) * cos(theta)) + (V * sin(alpha) * sin(theta));
zdot = -(V * cos(alpha) * sin(theta)) + (V * sin(alpha) * cos(theta));

Xdot = [V_dot,alpha_dot,q_dot,theta_dot,xdot,zdot]';
 
 end