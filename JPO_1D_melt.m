clc
clear;

%% Load geometry for plume equations
load('C:\Users\User\Desktop\matlab\Matlab\plume_1d_fem\realistic_geometry\section_1\realistic_geometry_sm300_section1.mat')
z = base_fl'; x = s_fl'; node = length(x); y=zeros(1,node);
Znot = z(1);


%% construct along-slope(flowline) coordinate y
 y(1) = 0; dx = zeros(1,node-1); dy=zeros(1,node-1);   dz =zeros(1,node-1);
 h=zeros(1,node-1);
for i=1:node-1
 %non-uniform coordinate
 dx(i) = x(i+1) - x(i);
 dz(i) = z(i+1) - z(i);
 
 dy(i) = sqrt( (dx(i))^2 + (dz(i))^2 )  ; 
 y(i+1) = y(i) + dy(i);
end

%% Define elments
h=dy;  

sin_theta = zeros(1,node-1); theta = zeros(1,node-1); deg_theta = zeros(1,node-1);
for i=1:node-1
 %% Minimal slope    
 sin_theta(i) = (z(i+1)-z(i))/(y(i+1)-y(i));  
 theta(i) = ( (asin(sin_theta(i)))); %radian
end


 %% temperature and salinity stratification (linear)
 temp_bottom = 0.47 ; % degree
 temp_grad = (-1.25 - temp_bottom)/600; % temperature decreaseing rate, degree/m 
 for i=1:node
 T_amb(i) = temp_grad*(z(i) - Znot) + temp_bottom;
 end
 sal_bottom = 34.73 ; %psu most saline
 sal_grad = (34 - sal_bottom)/600; %salinity decreaseing rate, psu/m 
 for i=1:node
 S_amb(i) = sal_grad*(z(i) - Znot) + sal_bottom;
 end

%% Variables for plume equations
 
Q = zeros(1,node); %flux, m^2/s
D = zeros(1,node); % thickness of plume, m
U = zeros(1,node); % plume velocity, m/s 
S = zeros(1,node); % plume salinity, psu
T = zeros(1,node); % plume temperature, C
QU = zeros(1,node);% momemtum flux, m3/s
QT = zeros(1,node); % temper flux, C*m2/s
QS = zeros(1,node); % salinity flux, psu*m2/s

%% Initial temperature and salinity of plume
lambda1 = -0.0573; %C/psu
lambda2 = 0.0832; %C
lambda3 = 7.61*10^-4; %C/m

% pressure freezing point at the groudning line depth with zero salinity (fresh)
S(1) = 0; 
T(1)  = lambda1*S(1) + lambda2 + lambda3*z(1);
 
%% Varying freshwater discharge
 width = 60; %comparable to the maximal diameter of the cone plume
 Q(1) = ((200*20)/width)*(3.6*10^(-2));%m^{2}s^{-1}
 QT(1) = Q(1)*T(1); QS(1) = Q(1)*S(1); 

%% parameters 
grav = 9.81;
E0 = 0.036; %entrainment coeffi, fixed for line plume
T_i =  -25;%C, temperature of the inner ice
beta_S = 7.86*10^-4; %1/psu
beta_T = 3.87*10^-5; %1/C
C_w = 3.974*10^3; %specific heat capacity of seawater, J/(kg*C)
C_i = 2.009*10^3; %specific heat capacity of ice, J/(kg*C)
L_i = 3.35*10^5; %latent heat of fusion for ice,  J/kg
CdTS0 = 6.0*10^-4;
Cd = 2.5*10^-3 ;%drag coeffi
trans_heat =2.2*10^-2; 
trans_sal = 6.2*10^-4; 
Cdtt = sqrt(Cd)*(trans_heat); % turbulent heat transfer coeffi
Cdss = sqrt(Cd)*(trans_sal); % turbulent salt transfer coeffi

%% Variables for the three meltrate equations 
T_bound = zeros(1,node); % boundary layer temperature
m_dot = zeros(1,node); % melting rate, m/s; positive for melting
e_dot =zeros(1,node);  %entrainment rate; m/s; e_dot = e_dot(U)
del_rho =zeros(1,node); % unitless; plume_density - ambient_water_density <0 for buoyancy 
 
del_rho(1) = beta_S*(S(1)-S_amb(1)) - beta_T*(T(1)-T_amb(1)); 

%% Initial conditions for plume velocity & thickness
% case1: setting U(1) assuming du/dx = 0 at x(1); balance velocity
numer = -Q(1)*del_rho(1)*grav*sin(theta(1));
denom = Cd + E0*sin(theta(1)); 
balance_velocity_line = (numer/denom)^(1/3) ; 

U(1) = balance_velocity_line;
D(1) =  Q(1)/U(1);
QU(1) = Q(1)*U(1);

%% y=1
[S_bound(1) T_bound(1) m_dot(1) myroot1 myroot2] = saline_bound(U(1),T(1),S(1),z(1)) ; 
e_dot(1) = E0*U(1)*sin(theta(1)); 
 
%% y=2 for U,Q,D,T,S
U(2) = U(1) ; %dU/dy = 0 at y=1 (Jenkins (1991))
Q(2) = Q(1) + dy(1)*(e_dot(1) + m_dot(1)) ; 
QU(2) = Q(2)*U(2);
D(2) = Q(2)/U(2);
  
RHS_QT =  ( e_dot(1)*T_amb(1) + m_dot(1)*T_bound(1) - Cdtt*U(1)*(T(1)-T_bound(1))   ) ;
QT(2) = QT(1) + dy(1)*RHS_QT;
T(2) = QT(2)/Q(2) ; 
  
RHS_QS = e_dot(1)*S_amb(1) + m_dot(1)*S_bound(1) - Cdss*U(1)*(S(1)-S_bound(1)) ;
QS(2) = QS(1) + dy(1)*RHS_QS;
S(2) = QS(2)/Q(2) ; 

%% Finite element method
K_L = [(-1/2 )  (1/2  ) ; (-1/2 )  (1/2  )];

%% Solve! (FEM + predictior, corrector)
for i =2:node-1

[S_bound(i) T_bound(i) m_dot(i) myroot11 myroot22] = saline_bound(U(i),T(i),S(i),z(i)) ; 

e_dot(i) = E0*U(i)*sin(theta(i)); 
del_rho(i) = beta_S*(S(i)-S_amb(i)) - beta_T*(T(i)-T_amb(i) ); 

%% Solve mass flux equ(FOR HAT)
RHS_Q_K1 = e_dot(i) + m_dot(i);
F_Q = [(RHS_Q_K1*(h(i)/2.) ); (RHS_Q_K1*(h(i)/2.) ) ];
Q_HAT = ( 0.5*F_Q(1,1) - K_L(1,1)*Q(i))/K_L(1,2);

%% Solve momentum flux equ(FOR HAT)
RHS_QU_K1 =  (   - D(i)*del_rho(i)*grav*sin(theta(i)) - Cd*(U(i)^2) );
F_QU = [(RHS_QU_K1*(h(i)/2.) ); (RHS_QU_K1*(h(i)/2.) ) ];
QU_HAT = ( 0.5*F_QU(1,1) - K_L(1,1)*QU(i))/K_L(1,2);

%% Update U&D(FOR HAT)
U_HAT = QU_HAT/Q_HAT;
D_HAT = Q_HAT/U_HAT;

%% Solve temperature flux equ & Update T(FOR HAT)
RHS_QT_K1 =  e_dot(i)*T_amb(i)   + m_dot(i)*T_bound(i) - Cdtt*abs(U(i))*(T(i)-T_bound(i)) ;
F_QT = [(RHS_QT_K1*(h(i)/2.) ); (RHS_QT_K1*(h(i)/2.) ) ];
QT_HAT = (0.5*F_QT(1,1) - K_L(1,1)*QT(i))/K_L(1,2);
T_HAT = QT_HAT/Q_HAT;

%% Solve salinity flux equ & Update S(FOR HAT)
RHS_QS_K1 =   e_dot(i)*S_amb(i) + m_dot(i)*S_bound(i) - Cdss*abs(U(i))*(S(i)-S_bound(i));
F_QS = [(RHS_QS_K1*(h(i)/2.) ); (RHS_QS_K1*(h(i)/2.) ) ];
QS_HAT = (0.5*F_QS(1,1) - K_L(1,1)*QS(i))/K_L(1,2);
S_HAT = QS_HAT/Q_HAT;


[S_bound_HAT T_bound_HAT m_dot_HAT myroot11 myroot22] = saline_bound(U_HAT,T_HAT,S_HAT,z(i)+dz(i)/2.0 ) ; 

e_dot_HAT = E0*U_HAT*sin(theta(i)); 
del_rho_HAT = beta_S*(S_HAT-S_amb(i)) - beta_T*(T_HAT-T_amb(i) ); 

%% Solve mass flux equ
RHS_Q_K2 = e_dot_HAT + m_dot_HAT;
F_Q_HAT = [(RHS_Q_K2*(h(i)/2.) ); (RHS_Q_K2*(h(i)/2.) ) ];
Q(i+1) = ( F_Q_HAT(1,1) - K_L(1,1)*Q(i))/K_L(1,2);

%% Solve momentum flux equ
RHS_QU_K2 =  ( - D_HAT*del_rho_HAT*grav*sin(theta(i)) - Cd*(U_HAT^2) );
F_QU_HAT = [(RHS_QU_K2*(h(i)/2.) ); (RHS_QU_K2*(h(i)/2.) ) ];
QU(i+1) = ( F_QU_HAT(1,1) - K_L(1,1)*QU(i))/K_L(1,2);

%% Update U&D
U(i+1) = QU(i+1)/Q(i+1);
D(i+1) = Q(i+1)/U(i+1);

%% Solve temperature flux equ
RHS_QT_K2 =  e_dot_HAT*T_amb(i) + m_dot_HAT*T_bound_HAT - Cdtt*abs(U_HAT)*(T_HAT -T_bound_HAT);
F_QT_HAT = [(RHS_QT_K2*(h(i)/2.) ); (RHS_QT_K2*(h(i)/2.) ) ];
QT(i+1) = ( F_QT_HAT(1,1) - K_L(1,1)*QT(i))/K_L(1,2);
T(i+1) = QT(i+1)/Q(i+1);

%% Solve salinity flux equ
RHS_QS_K2 =   e_dot_HAT*S_amb(i) + m_dot_HAT*S_bound_HAT - Cdss*abs(U_HAT)*(S_HAT-S_bound_HAT);%hat
F_QS_HAT = [(RHS_QS_K2*(h(i)/2.) ); (RHS_QS_K2*(h(i)/2.) ) ];
QS(i+1) = ( F_QS_HAT(1,1) - K_L(1,1)*QS(i))/K_L(1,2);
S(i+1) = QS(i+1)/Q(i+1); 



%% Stop for integration (either negative buoyancy | neutral buoyancy)
if  ( sin_theta(i-1)<=0 ) | del_rho(i-1)*del_rho(i)<0
break
end

end
run1 = i;

meltrate_per_year = (365*24*3600).*m_dot;
mean_meltrate_per_year = mean(meltrate_per_year);
maximum_meltrate_per_year = max(meltrate_per_year);
mininum_meltrate_per_year = min(meltrate_per_year);

%% coarser coordinate for parameterization (don't have to use finer resolutions as used for plume model)
  
load('C:\Users\User\Desktop\matlab\Matlab\plume_1d_fem\realistic_geometry\section_1\realistic_geometry_sm300_coarser_section1.mat')
zc = base_fl'; xc = s_fl'; nodec = length(xc); yc=zeros(1,nodec);
Znot = zc(1);

%% construct along-slope(flowline) coordinate Y
 yc(1) = 0; dxc = zeros(1,nodec-1); dyc=zeros(1,nodec-1);   dzc =zeros(1,nodec-1);
 hc=zeros(1,nodec-1);

 for i=1:nodec-1
 dxc(i) = xc(i+1) - xc(i);
 dzc(i) = zc(i+1) - zc(i);
 dyc(i) = sqrt( (dxc(i))^2 + (dzc(i))^2 )  ; 
 yc(i+1) = yc(i) + dyc(i);
end

%% Define elments
hc=dyc; %elements 

sin_thetac = zeros(1,nodec-1); thetac = zeros(1,nodec-1); 
for i=1:nodec-1
 %% Minimal slope    
 sin_thetac(i) = (zc(i+1)-zc(i))/(yc(i+1)-yc(i));  
 thetac(i) = ( (asin(sin_thetac(i)))); %radian
 tan_thetac(i) = (zc(i+1)-zc(i))/(xc(i+1)-xc(i));
 %sin_thetac2(i) = sin(atan(tan_thetac(i)));

end
 
%% temperature and salinity stratification 
 temp_bottom = 0.47 ; %psu
 temp_grad = (-1.25 - temp_bottom)/600;%; 
 for i=1:nodec
 T_ambc(i) = temp_grad*(zc(i) - Znot) + temp_bottom;
 end
 sal_bottom = 34.73 ; %psu most saline
 sal_grad = (34 - sal_bottom)/600; 
 for i=1:nodec
 S_ambc(i) = sal_grad*(zc(i) - Znot) + sal_bottom;
 end

%% Initial values 

Sc(1) = 0;
Tc(1)  = lambda1*Sc(1) + lambda2 + lambda3*zc(1);
del_rhoc =zeros(1,nodec); % unitless; plume_density - ambient_water_density <0 for buoyancy 
del_rhoc(1) = beta_S*(Sc(1)-S_ambc(1)) - beta_T*(Tc(1)-T_ambc(1)); 
T_freezing = lambda1*Sc + lambda2 + lambda3*zc;
T_m = (C_i/C_w)*T_i - L_i/C_w ;

%% Source correction parameter 
ystar = ((Q(1)^2)/((E0^2)*(  (  sin_theta(1)  )^3)*(-del_rhoc(1)*grav)))^(1/3);
zstar = ystar*sin_theta(1) ;
rho_sw = 1028; % density of seawater kg/m^3
corrected_depth = Znot - zstar;

%% Finite source correction
% T_melting_corr is the pressure freezing point of freshwater at the depth of the virtual source
% Note that the pressure exerted should be calculated using the density of seawater
% evaluated using the Clausius- Clapeyron equation
% lambda 2 + lambda3*corrected_depth is the pressure freezing point of seawater at the depth of the virtual source
% evaluated using the linearized Clausius- Clapeyron equation 

pressure_corr = (rho_sw)*grav*( - corrected_depth );
T_melting_corr =  ( 273.15*exp( -0.00027196*(pressure_corr*10^-6 - 0.1) ) ) - 273.15  % 0.1 is for atmospheric pressure in the unit of MPa

% virtual thermal forcing \Delta T_{z = Z_{gl} + zstar}
TF_source = T_melting_corr - (lambda2 + lambda3*( corrected_depth )) ; 


%% parameter for unitless thermal forcing 
aa = -0.6635;
bb =  0.2581;
cc1 =  1;
cc2 = 0.0065;
cc3 = -0.0022;
rr = 0.53119;
dd = -0.0486;
ee = 0.0651;
A1 = 2.1449 ; A2  = -0.0094 ;  b1 = -0.1003;
alpha_u0 =  0.5585;
alpha_u1 = -0.6858;
mean_angle = mean(sin_thetac);
GL_angle = max(sin_thetac);

%% parameterization calculation
for i=1:nodec-1

   T_af(i) = lambda1*S_ambc(i) + lambda2 + lambda3*zc(i); % The ambient freezing point 
   M0_app(i) = CdTS0./(T_af(i) - T_m); % approximation of melt rate factor, Eq. (14)
   balance_thermal_forcing(i) = ( (T_ambc(i) - T_af(i)) )...
                                /(1 - ((M0_app(i)*T_m)/( E0*sin_thetac(i)  ))); % Balance thermal forcing, Eq. (16)
   balance_thermal_forcing_max_angle(i) = ( (T_ambc(i) - T_af(i)) )...
                               /(1 -  ((M0_app(i)*T_m)/( E0*GL_angle   ))); %the balance thermal forcing that is defined using the grounding line slope angle:
   %% geometrical temperature factor-weighted average 
   epsilon_mean = (E0*mean_angle)/(CdTS0 + E0*mean_angle);
   epsilon_gl = (E0*GL_angle)/(CdTS0 + E0*GL_angle);
   epsilon_coffi = 2; 
  
   %% Formulation 43
   balance_thermal_forcing_weighted(i) =   ( (epsilon_gl^epsilon_coffi).*balance_thermal_forcing(i) + (epsilon_mean^epsilon_coffi).*balance_thermal_forcing_max_angle(i) )./(epsilon_gl^epsilon_coffi + epsilon_mean^epsilon_coffi );
   ambient_thermal_forcing_samb(i) = T_ambc(i) - T_af(i);% ambient thermal forcing  

%% Unitless thermal forcing parameter
GT(i) =   ( 1 )/( 1 -  (M0_app(i)*T_m)/(E0*sin_thetac(i))  );
alpha_T = 1.2182*(ambient_thermal_forcing_samb(1)^-0.0231)*(GT(1)^-0.72);
KM = ((1*Q(1))^(cc2*balance_thermal_forcing(1)^rr  + cc3) )* ( (1*balance_thermal_forcing(1))...
                        ^(0.67*(Q(1))^0.2553 + aa*Q(1)^bb +dd/balance_thermal_forcing(1) + ee ) );

%% Test stratification using the refined form where changing basal slope angle is accounted for 
factor_param = KM*A1*(abs(Znot)^b1)*( 1. + A2*ambient_thermal_forcing_samb(1) )*(M0_app(1));
balance_thermal_forcing_final(i) = factor_param*( balance_thermal_forcing_weighted(i)*( ...
                                  (1 - ((( ystar )./(yc(i) + ystar  -yc(1)  )  )^(alpha_T)  ) )  ) ) ;
%Eq. 35
source_correction(i) = factor_param*( TF_source*( ...
                                  (1 - ((( ystar )./(yc(i) + ystar  -yc(1)  )  )^(alpha_T)  ) )  ) ) ;

%Eq. 46
balance_thermal_forcing_final(i) =  balance_thermal_forcing_final(i) + source_correction(i);

%% Refined velocity for changing slope (eq. 39-42)
alpha_u(i) = alpha_u0.*( sin_thetac(i) ).^(alpha_u1);  %same form for constant slope
numer(i) = -Q(1)*del_rhoc(1)*grav*sin_thetac(i);
denom(i) = Cd + E0*sin_thetac(i); 
b_m(i) = grav*sin_thetac(i)*(beta_S*S_ambc(i) - beta_T*(T_ambc(i) - T_m ));
balance_velocity(i) = sqrt(  (numer(i)/denom(i))^(2/3) ...
                        + M0_app(i)*(2/3)*(b_m(i)/(Cd + E0*sin_thetac(i)  ))*( ...
                          balance_thermal_forcing(1)*yc(i)  ));
balance_velocity_final(i) = U(1) + (balance_velocity(i) - U(1))*(1 -  (ystar./(yc(i) + ystar))^(alpha_u(i))   );
%% parameterized melt rates for changing slope (eq. 48 )
para_meltrate_final(i) = balance_velocity_final(i)*balance_thermal_forcing_final(i);
end 



figure(1)
plot(y(1:run1)./10^3,365*24*3600.*m_dot(1:run1),'black-','linewidth',1)
hold on
plot(yc(1:length(para_meltrate_final))./10^3,365*24*3600.*para_meltrate_final,'red-','linewidth',1)
xlim([0 y(end)/1000])
legend('Plume model','Parameterization')
xlabel('x (km)')
ylabel('$\dot m$ (m/yr)','Interpreter','latex')


figure(2)
plot(y(1:run1)./10^3,U(1:run1),'black-','linewidth',1)
hold on
plot(yc(1:length(balance_velocity_final))./10^3,balance_velocity_final,'red-','linewidth',1)
xlim([0 y(end)/1000])
legend('Plume model','Parameterization')
xlabel('x (km)')
ylabel('$U$ (m/s)','Interpreter','latex')





