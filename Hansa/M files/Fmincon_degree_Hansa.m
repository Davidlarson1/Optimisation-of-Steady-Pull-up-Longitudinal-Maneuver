clear all;
clc;
format short;


 options = optimoptions('fmincon','StepTolerance',1e-4, 'MaxFunctionEvaluations', 10000,'ConstraintTolerance',1e-10,'MaxFunEvals',10000,'MaxIter',50000);
 [q, fval, exitflag, output] = fmincon(@objfun,[0.8],[],[],[],[],[0],[2],@mycons,options); 

 global Ve t  N q X T de n a_1 M n_1 alpha_trim T_trim AeroPara GeoPara de_trim gamma TF TF_trim

  
% AOA
subplot(5,1,1);
plot(t,X(1,:) ,'k-','LineWidth',1);
xlabel('t(s)');
ylabel('\alpha (deg)');
title('AOA Vs t');
grid on

% Velocity
%subplot(4,1,1);
%plot(t,Ve* ones(1,N),'k-','LineWidth',1);
%xlabel('t(s)');
%ylabel('V (m/s)');
%title('Velocity Vs t');
grid on

%Pitch rate q
%subplot(4,1,3);
%plot(t,(q * ones(1,N)) ,'k-','LineWidth',1);
%xlabel('t(s)');
%ylabel('q(deg/s)');
%title('Pitch Rate Vs t');
grid on

% Load Factor
subplot(5,1,3);
plot(t,n,'k-','LineWidth',1);  
xlabel('t(s)');
ylabel('n');
title('Load Factor Vs t');
grid on

%Thrust Factor
subplot(5,1,2);
plot(t,TF,'k-','LineWidth',1);
xlabel('t(s)');
ylabel('TF');
title(' Thrust Factor Vs t');
grid on

% Elevator
subplot(5,1,4);
plot(t,de ,'k-','LineWidth',1);%* 180/pi
xlabel('t(s)');
ylabel('\delta_e(deg)');
title('Elevator Vs t');
grid on

% Body Z
%subplot(4,1,4);%
%plot(X(2,:),-X(3,:),'k-','LineWidth',1);
%xlabel('Time(s)');
%xlabel('X(m)');
%ylabel('Altitude(m)');
%title(' X Vs Altitude');
grid on

% Flight path Angle
subplot(5,1,5);
plot(t,gamma ,'k-','LineWidth',1);
xlabel('t(s)');
ylabel('\gamma (deg)');
title('Flight Path Angle Vs t');
grid on
   
        
% Objective Function

function [f] = objfun(q)
      
        f = -q^2;
end

function [c, ceq] = mycons(q)

global Ve t  N  X  n T de  M a_1 n_1 de_trim alpha_trim T_trim AeroPara GeoPara gamma TF TF_trim

  
     Ve = 40;     %input('Enter desired Velocity  :   ');
     TF_max =9.15;% input('Enter maximum Thrust  :   ');
     alpha_max = 66 ;%input('Enter Stall Alpha  :   ');
     n_max = 4;% input('Enter maximum load factor  :   ');
     
AeroPara = [0.354,0.036,0.052,4.972 * 0.0175 , 37.259 * 0.0175 ,-1.008 * 0.0175, -0.496 * 0.0175,-11.729 * 0.0175,0.265 * 0.0175 ,0.026 * 0.0175,0.061* 0.0175 ];
 GeoPara = [750,9.81 ,1.225,12.47, 1.21, 0.16];

    
    mass = GeoPara(1);g = GeoPara(2);rho = GeoPara(3);S = GeoPara(4);c_bar = GeoPara(5);k = GeoPara(6);
    CLO = AeroPara(1);CDO = AeroPara(2);CmO = AeroPara(3);
    CL_alpha = AeroPara(4); CL_q = AeroPara(5);
    Cm_de = AeroPara(6); Cm_alpha = AeroPara(7);Cm_q = AeroPara(8);CL_de = AeroPara(9);CD_de = AeroPara(10);CD_alpha = AeroPara(11);

     
     W = mass * g;

CL_trim = (2 * W)/(rho * S * Ve^2);
CD_trim = CDO + ( CD_alpha * CL_trim^2);

alpha_trim = (CL_de*CmO - CLO*Cm_de + CL_trim*Cm_de)/(CL_alpha*Cm_de - CL_de*Cm_alpha);
de_trim = -(CL_alpha*CmO - CLO*Cm_alpha + CL_trim*Cm_alpha)/(CL_alpha*Cm_de - CL_de*Cm_alpha);
T_trim  = W/(CL_trim/CD_trim);
TF_trim = T_trim/W;

     % Time changes
     
     dt = 0.1;
     tf = 180/q;
     t = 0:dt:tf;
     N = (length(t));
     
  
    x0 = [alpha_trim,0,-500,alpha_trim]';
    X = zeros(4, N);
    X(:,1) = x0;
    
    de = zeros(1, N);
    de(:,1) = de_trim;
    n = zeros(1,N);
    T = zeros(1,N);
    T(:,1) = T_trim;
    n(:,1) = 0;
     TF = zeros(1,N);
    TF(:,1) = TF_trim;
    gamma = zeros(1,N); % Flight path angle = theta-alpha
    gamma(:,1) = 0;
    

    a = [q,Ve]';
    
  
    
   for i = 1 : N-1
       
           
        k1 =  dt * Equations_Degree((X(:,i)),a);
        
        k2 = dt * Equations_Degree((X(:,i)  + (.5 * k1)),a);
          
        k3 = dt *  Equations_Degree((X(:,i) + (.5 * k2)),a);
            
        k4 = dt *  Equations_Degree((X(:,i) +(k3)),a);
          
        X(:,i+1) = X(:,i) + ((k1 + 2*k2 + 2*k3 + k4)/6);
        
        gamma(1,i +1) = X(4,i +1) - X(1,i +1);
        
        CL(1,i+1) = CLO + (CL_alpha * X(1,i+1)) + (CL_q *(q * c_bar)/(2*Ve));
        CD(1,i+1) = (CDO + (CD_alpha * X(1,i+1)));
        
        L(1,i+1) = .5 * rho * Ve^2 * S * CL(1,i+1);
        D(1,i+1) = .5 * rho * Ve^2 * S * CD(1,i+1);
        
        n(1,i+1) = L(1,i+1)/(mass * g);
        de(1,i +1) = -(CmO + (Cm_alpha * X(1,i+1)) + (Cm_q *((q * c_bar)/(2*Ve))))/(Cm_de);
        
        % This gives the value of Thrust or engine power
        T(1,i +1)  = (D(1,i+1) + (mass * g * sind(X(4,i+1) - X(1,i+1))))/cosd(X(1,i+1));
        TF(1,i+1) = T(1,i+1)/W;

   end
   
% Limit for Thrust/AOA/LOAD factor
     n_1 =  max(n, [], 'all');
      M =  max(TF, [], 'all'); 
      a_1 =  max(X(1,:), [], 'all'); 
     % de_1 = min(de(1,:), [], 'all'); 
     % de_2 = max(de(1,:), [], 'all');

      c = [a_1 - alpha_max, M - TF_max, n_1 - n_max]; % [a_1 - alpha_max, M - T_max, n_1 - n_max];
      ceq = [];
  
end
