clear all;
clc;
format short;


 options = optimoptions('fmincon','StepTolerance',1e-6, 'MaxFunctionEvaluations', 1000,'ConstraintTolerance',1e-15,'MaxFunEvals',10000,'MaxIter',50000);
 [q, fval, exitflag, output] = fmincon(@objfun,[0.1],[],[],[],[],[0],[3],@mycons,options); 

 global Ve t  N  X T n TF LF de a_1 M n_1 alpha_trim AeroPara GeoPara de_trim gamma
  
      
% AOA
subplot(5,1,1);
plot(t,X(1,:)* 180/pi ,'k-','LineWidth',1);
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
%plot(t,(q * ones(1,N)) * 180/pi ,'k-','LineWidth',1);
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
plot(t,de * 180/pi ,'k-','LineWidth',1);%* 180/pi
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
plot(t,gamma * 180/pi,'k-','LineWidth',1);
xlabel('t(s)');
ylabel('\gamma (deg)');
title('Flight Path Angle Vs t');
grid on
  
   
        
% Objective Function

function [f] = objfun(q)
        % Objective function
        f = -q^2;
end

function [c, ceq] = mycons(q)

global Ve t  N  X T de a_1 M n_1 alpha_trim AeroPara GeoPara gamma de_trim n TF


     Ve = 150;     %input('Enter desired Velocity  :   ');
    
     % Time changes
     
     dt = 0.1;
     tf = pi/q;
     t = 0:dt:tf;
     N = (length(t));
     
   
 GeoPara = [1128.09*14.593 , 9.81 , 1.225 , 400*(0.305^2) , 11.52*0.305 , 239720.815 , 239720.815 ,259969.9570048 ,-2.6, 37.42*0.305];
  mass = GeoPara(1);g = GeoPara(2);rho = GeoPara(3);S = GeoPara(4);c_bar = GeoPara(5);

%   Aero parameters for a constant alpha = 2.097 deg 
  AeroPara0 = Aero(0.0349);
 % AeroPara0 = [0.0108667726033608,-0.166894131707912,-4.88638129287557,-0.855509415408481,-0.0481837305138529,5.81591056833472,4.37262151305434,0.775933175385267,0.00913882467000680,0.498377553216749,-0.0466353533722922,-0.125567837094329];

 CL0 = AeroPara0(5);CD0= AeroPara0(9);Cm0 = AeroPara0(1);
Cm_alpha = AeroPara0(2);CL_alpha = AeroPara0(6);CD_alpha = AeroPara0(10);
CL_q= AeroPara0(7) ;CD_q = AeroPara0(11);Cm_q = AeroPara0(3);
Cm_del = AeroPara0(4); CL_del = AeroPara0(8);CD_del = AeroPara0(12);

%AeroPara0 = [Cm0,Cm_alpha,Cm_q,Cm_del,CL0,CL_alpha ,CL_q,CL_del,CD0,CD_alpha,CD_q,CD_del];


 % Trim Calculation
 
   % Trim Calculation
W = mass * g;

CL_trim = (2 * W)/(rho * S * Ve^2);

alpha_trim = (CL_del*Cm0 - CL0*Cm_del + CL_trim*Cm_del)/(CL_alpha*Cm_del - CL_del*Cm_alpha);
de_trim = -(CL_alpha*Cm0 - CL0*Cm_alpha + CL_trim*Cm_alpha)/(CL_alpha*Cm_del - CL_del*Cm_alpha);

CD_trim =  CD0 + (CD_alpha * alpha_trim) ; 
T_trim  = W/(CL_trim/CD_trim);
L_trim = 0.5 * rho * Ve^2 * S * CL_trim;
D_trim = 0.5 * rho * Ve^2 * S * CD_trim;

    AeroPara0 = Aero(alpha_trim );
    
    x0 = [alpha_trim ,0,-500,alpha_trim ]';
    X = zeros(4, N);
    n = zeros(1,N);
    n(:,1) = L_trim/W;
    X(:,1) = x0;
    
    AeroPara = zeros(12, N);
    AeroPara(:,1) = AeroPara0;
    
    de = zeros(1, N);
    de(:,1) = de_trim ;
    T = zeros(1,N);
    TF = zeros(1,N);
    TF(:,1) = T_trim/W;
    T(:,1) = T_trim;
   
    CL = zeros(1,N);
    CD = zeros(1,N);
    L = zeros(1,N);
    D = zeros(1,N);
    CL(:,1) = CL_trim;
    L(:,1) = L_trim;
    CD(:,1) = CD_trim;
    D(:,1) = D_trim;
    
    gamma = zeros(1,N);
    gamma(:,1) = 0;

   a = [q,Ve]';
    
 
   for i = 1 : N-1
       
           
        k1 =  dt * Equations_F18((X(:,i)),a,AeroPara(:,i));
        
        k2 = dt * Equations_F18((X(:,i)  + (.5 * k1)),a,AeroPara(:,i));
          
        k3 = dt *  Equations_F18((X(:,i) + (.5 * k2)),a,AeroPara(:,i));
            
        k4 = dt *  Equations_F18((X(:,i) +(k3)),a,AeroPara(:,i));
          
        X(:,i+1) = X(:,i) + ((k1 + 2*k2 + 2*k3 + k4)/6);
        
        gamma(1,i+1) = X(4,i+1) - X(1,i+1);
        
        AeroPara(:,i+1)= Aero(X(1,i+1));

        de(1,i +1) = -(AeroPara(1,i+1) + AeroPara(2,i+1)*X(1,i+1) +(AeroPara(3,i+1)*((q * c_bar)/(2*Ve))))/(AeroPara(4,i+1));
        
        CL(1,i+1) = AeroPara(5,i+1) + AeroPara(6,i+1)*X(1,i+1) + (AeroPara(7,i+1) *(q * c_bar)/(2*Ve)) + AeroPara(8,i+1) * de(1,i+1);
        CD(1,i+1) = (AeroPara(9,i+1) + AeroPara(10,i+1)*X(1,i+1) + AeroPara(12,i+1) * de(1,i+1));
        
        L(1,i+1) = .5 * rho * Ve^2 * S * CL(1,i+1);
        D(1,i+1) = .5 * rho * Ve^2 * S * CD(1,i+1);
        
        n(1,i+1) = L(1,i+1)./(W);
     
        T(1,i+1)  = (D(1,i+1) + (mass * g * sin(X(4,i+1) - X(1,i+1))) )/cos(X(1,i+1));
        TF(1,i+1) =  T(1,i+1)/(W);

   end
   
 %Limit for Thrust/AOA/LOAD factor
     n_1 =  max(n, [], 'all');
      M =  max(TF, [], 'all'); 
      a_1 = 57.29 *  max(X(1,:), [], 'all'); 
     % de_1 = min(de(1,:), [], 'all'); 
     % de_2 = max(de(1,:), [], 'all');

      c = [a_1  - 40, M - 1.18, n_1 - 6]; 
      ceq = [];
    

end