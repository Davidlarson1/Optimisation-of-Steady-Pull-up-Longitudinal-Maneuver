
function X_dot = Equations_F18(X0,a,AeroPara0)

alpha_1 = X0(1); 
x = X0(2);
z = X0(3);
theta_1 = X0(4); 

q_1 = a(1);Ve_1 = a(2);


   GeoPara = [1128.09*14.593 , 9.81 , 1.225 , 400*(0.305^2) , 11.52*0.305 , 239720.815 , 239720.815 ,259969.9570048 ,-2.6, 37.42*0.305];
   mass = GeoPara(1);g = GeoPara(2);rho = GeoPara(3);S = GeoPara(4);c_bar = GeoPara(5);

   
   

CL0 = AeroPara0(5);CD0= AeroPara0(9);Cm0 = AeroPara0(1);
Cm_alpha = AeroPara0(2);CL_alpha = AeroPara0(6);CD_alpha = AeroPara0(10);
CL_q= AeroPara0(7) ;CD_q = AeroPara0(11);Cm_q = AeroPara0(3);
Cm_del = AeroPara0(4); CL_del = AeroPara0(8);CD_del = AeroPara0(12);

AeroPara0 = [Cm0,Cm_alpha,Cm_q,Cm_del,CL0,CL_alpha ,CL_q,CL_del,CD0,CD_alpha,CD_q,CD_del];


% Aerodynamic and Geometric Parameters

CL1 = (CL0 + (CL_alpha * alpha_1)+ (CL_q *((q_1 * c_bar)/(2*Ve_1))));
CD1 = (CD0 + (CD_alpha * alpha_1)+ (CD_q *((q_1 * c_bar)/(2*Ve_1)))); 
Cm1 = (Cm0 + (Cm_alpha * alpha_1)+ (Cm_q *((q_1 * c_bar)/(2*Ve_1)))); 


alpha_dot = -(rho * Ve_1* S)/(2 * mass) * (CL1 + (CD1 * tan(alpha_1)) - (Cm1 *  (CL_del +  (CD_del * tan(alpha_1)))/ (Cm_del ))) + (g * cos(theta_1))/(Ve_1 * cos(alpha_1)) + q_1;
xdot = (Ve_1 * cos(alpha_1) * cos(theta_1)) + (Ve_1 * sin(alpha_1) * sin(theta_1));
zdot = -(Ve_1 * cos(alpha_1) * sin(theta_1)) + (Ve_1 * sin(alpha_1) * cos(theta_1));
theta_dot = q_1;

X_dot = [alpha_dot,xdot,zdot,theta_dot]';

 end