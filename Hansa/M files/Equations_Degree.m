
function X_dot = Equations_Degree(X0,a)
 
alpha_1 = X0(1); 
x = X0(2);
z = X0(3);
theta_1 = X0(4); 

q_1 = a(1);Ve_1 = a(2);

GeoPara = [750,9.81 ,1.225,12.47, 1.21, 0.16];
mass = GeoPara(1);g = GeoPara(2);rho = GeoPara(3);S = GeoPara(4);c_bar = GeoPara(5);k = GeoPara(6);


AeroPara = [0.354,0.036,0.052,4.972 * 0.0175 , 37.259 * 0.0175 ,-1.008 * 0.0175, -0.496 * 0.0175,-11.729 * 0.0175,0.265 * 0.0175,0.026 * 0.0175,0.061* 0.0175 ];

CLO = AeroPara(1);CDO = AeroPara(2);CmO = AeroPara(3);
CL_alpha = AeroPara(4); CL_q = AeroPara(5);
Cm_de = AeroPara(6); Cm_alpha = AeroPara(7);Cm_q = AeroPara(8);CL_de = AeroPara(9);CD_de = AeroPara(10);CD_alpha = AeroPara(11);

CL1 = (CLO + (CL_alpha * alpha_1) + (CL_q *((q_1 * c_bar)/(2*Ve_1))));
CD1 = (CDO + (CD_alpha * alpha_1));
Cm1 = (CmO + (Cm_alpha * alpha_1) + (Cm_q *((q_1 * c_bar)/(2*Ve_1))));

alpha_dot = -(rho * Ve_1* S)/(2 * mass) * [CL1 + (CD1 * tand(alpha_1)) - (Cm1 * (CL_de + CD_de * tand(alpha_1))/(Cm_de))] + (g * cosd(theta_1))/(Ve_1 * cosd(alpha_1)) + q_1;
xdot = (Ve_1 * cosd(alpha_1) * cosd(theta_1)) + (Ve_1 * sind(alpha_1) * sind(theta_1));
zdot = -(Ve_1 * cosd(alpha_1) * sind(theta_1)) + (Ve_1 * sind(alpha_1) * cosd(theta_1));
theta_dot = q_1;

X_dot = [alpha_dot,xdot,zdot,theta_dot]';

 end