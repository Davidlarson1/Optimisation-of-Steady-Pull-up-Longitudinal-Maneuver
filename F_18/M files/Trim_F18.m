
clear all;
clc;
format short;

  Ve = 220; % Trim velocity
    
   GeoPara = [1128.09*14.593 , 9.81 , 1.225 , 400*(0.305^2) , 11.52*0.305 , 239720.815 , 239720.815 ,259969.9570048 ,-2.6, 37.42*0.305];
   mass = GeoPara(1);g = GeoPara(2);rho = GeoPara(3);S = GeoPara(4);c_bar = GeoPara(5);

   AeroPara0 = Aero(0.0349);
   
 CL0 = AeroPara0(5);CD0= AeroPara0(9);Cm0 = AeroPara0(1);
Cm_alpha = AeroPara0(2);CL_alpha = AeroPara0(6);CD_alpha = AeroPara0(10);
CL_q= AeroPara0(7) ;CD_q = AeroPara0(11);Cm_q = AeroPara0(3);
Cm_del = AeroPara0(4); CL_del = AeroPara0(8);CD_del = AeroPara0(12);

 % Trim Calculation
W = mass * g;

CL_trim = (2 * W)/(rho * S * Ve^2);

alpha_trim = (CL_del*Cm0 - CL0*Cm_del + CL_trim*Cm_del)/(CL_alpha*Cm_del - CL_del*Cm_alpha);
de_trim = -(CL_alpha*Cm0 - CL0*Cm_alpha + CL_trim*Cm_alpha)/(CL_alpha*Cm_del - CL_del*Cm_alpha);

CD_trim =  CD0 + (CD_alpha * alpha_trim) + CD_del * de_trim ; 
T_trim  = W/(CL_trim/CD_trim);

q0 = 0;

    % Time
    
     dt = 0.1;
     tf = 1200;
     t = 0:dt:tf;
     N = (length(t));
     
    de =  de_trim ;
    T = T_trim ;
    
    x0 = [Ve,alpha_trim,q0,alpha_trim,0,-500]';
    X = zeros(6, N);
    X(:,1) = x0;
  
    
   for i = 1 : N-1
       
           
        k1 =  dt * Trim_equations_F18((X(:,i)),T,de,AeroPara0);
        
        k2 = dt * Trim_equations_F18((X(:,i)  + (.5 * k1)),T,de,AeroPara0);
          
        k3 = dt *  Trim_equations_F18((X(:,i) + (.5 * k2)),T,de,AeroPara0);
            
        k4 = dt *  Trim_equations_F18((X(:,i) +(k3)),T,de,AeroPara0);
          
        X(:,i+1) = X(:,i) + ((k1 + 2*k2 + 2*k3 + k4)/6);
        
        q_Round = round(X(3,:),6); % Rounding the value of pitch rate q.
    
        theta_Round = round(X(4,:),4); % Rounding the value of pitch rate theta.
        

   end
   
% AOA
subplot(4,1,2);
plot(t,X(2,:)* 180/pi,'k-','LineWidth',1);
xlabel('time(s)');
ylabel('\alpha(deg)');
title('AOA Vs t');
grid on

% Velocity
subplot(4,1,1);
plot(t,X(1,:),'k-','LineWidth',1);
xlabel('time(s)');
ylabel('V(m/s)');
title('Velocity Vs t');
grid on

% Pitch rate q
subplot(4,1,3);
plot(t,q_Round* 180/pi,'k-','LineWidth',1);
xlabel('time(s)');
ylabel('q(deg/s)');
title('Pitch rate Vs t');
grid on

%Theta
subplot(4,1,4);
plot(t,X(4,:)* 180/pi,'k-','LineWidth',1);  
xlabel('time(s)');
ylabel('\theta(deg)');
title('Pitch Angle Vs t');
grid on

