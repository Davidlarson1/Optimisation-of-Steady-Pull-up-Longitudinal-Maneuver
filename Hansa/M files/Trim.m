
clear all;
clc;
format short;

  Ve = 40; % Trim velocity
    
  GeoPara = [750,9.81 ,1.225,12.47, 1.21, 0.16];
  mass = GeoPara(1);g = GeoPara(2);rho = GeoPara(3);S = GeoPara(4);c_bar = GeoPara(5);

AeroPara = [0.354,0.036,0.052,4.972 * 0.0175 , 37.259 * 0.0175 ,-1.008 * 0.0175, -0.496 * 0.0175,-11.729 * 0.0175,0.265 * 0.0175,0.026 * 0.0175,0.061* 0.0175 ];

CL0 = AeroPara(1);CD0 = AeroPara(2);Cm0 = AeroPara(3);
CL_alpha = AeroPara(4); CL_q = AeroPara(5);
Cm_de = AeroPara(6); Cm_alpha = AeroPara(7);Cm_q = AeroPara(8);CL_de = AeroPara(9);CD_de = AeroPara(10);CD_alpha = AeroPara(11);

 % Trim Calculation
W = mass * g;

CL_trim = (2 * W)/(rho * S * Ve^2);

alpha_trim = (CL_de*Cm0 - CL0*Cm_de + CL_trim*Cm_de)/(CL_alpha*Cm_de - CL_de*Cm_alpha);
de_trim = -(CL_alpha*Cm0 - CL0*Cm_alpha + CL_trim*Cm_alpha)/(CL_alpha*Cm_de - CL_de*Cm_alpha);

CD_trim =  CD0 + (CD_alpha * alpha_trim) + CD_de * de_trim ; 
T_trim  = W/(CL_trim/CD_trim);
q0 = 0;

    % Time
    
     dt = 0.1;
     tf = 400;
     t = 0:dt:tf;
     N = (length(t));
     
    de =  de_trim ;
    T = T_trim ;
    
    x0 = [Ve,alpha_trim,0,alpha_trim,0,-500]';
    X = zeros(6, N);
    X(:,1) = x0;
  
    
   for i = 1 : N-1
       
           
        k1 =  dt * Trim_equation((X(:,i)),T,de);
        
        k2 = dt * Trim_equation((X(:,i)  + (.5 * k1)),T,de);
          
        k3 = dt *  Trim_equation((X(:,i) + (.5 * k2)),T,de);
            
        k4 = dt *  Trim_equation((X(:,i) +(k3)),T,de);
          
        X(:,i+1) = X(:,i) + ((k1 + 2*k2 + 2*k3 + k4)/6);
        
        q_Round = round(X(3,:),6); % Rounding the value of pitch rate q.
    
        theta_Round = round(X(4,:),4); % Rounding the value of pitch rate theta.
        

   end
   
% AOA
subplot(4,1,2);
plot(t,X(2,:),'k-','LineWidth',1);
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
plot(t,q_Round,'k-','LineWidth',1);
xlabel('time(s)');
ylabel('q(deg/s)');
title('Pitch rate Vs t');
grid on

%Theta
subplot(4,1,4);
plot(t,X(4,:),'k-','LineWidth',1);  
xlabel('time(s)');
ylabel('\theta(deg)');
title('Pitch Angle Vs t');
grid on

% Thrust
%subplot(4,1,2);
%plot(t,T * ones(1,N),'k-','LineWidth',1);
%xlabel('time(s)');
%ylabel('Thrust(N)');
%title('Thrust Vs t');
grid on

% de_trim
%subplot(4,1,1);
%plot(t,de * ones(1,N),'k-','LineWidth',1);
%xlabel('time(s)');
%ylabel('\delta_e(deg)');
%title('\delta_e Vs t');
grid on

% XZ
%subplot(4,1,3);
%plot(X(5,:),X(6,:),'k-','LineWidth',1);
%xlabel('X(m)');
%ylabel('Z(m)');
%title('X Vs Z');
grid on

%Dummy
%subplot(4,1,4);
%plot(t,X(4,:),'k-','LineWidth',1);  
%xlabel('time(s)');
%ylabel('\theta(deg)');
%title('\theta Vs t');
grid on