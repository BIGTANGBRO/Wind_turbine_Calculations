clc
clear

% This script considers the following the structure strength of airfoil
% sg6043;
[Cl,Cd,c,beta,r,step,R]=GetData();
beta_angle=beta;
beta=beta*pi/180;%turn the beta angle to the radian
% Input Air
Vwind=12; % U infinity - Wind Speed; 6,8,10,12 ms^-1
omega_rpm=input('input the angular velocity of the turbine'); % Angular velocity; maximum 3000 rpm
a1=1/3; % Induction factor
rho=1.225; % Density; kg/m^3

% Input Geometry
nb=input('input the number of the blades:'); % Number of blades
np=length(c); % Number of panels
% Material Properties
E=2230*10^6; % Tensile modulus; Pa
eps_max=0.07; % Break elongation
sigma_max=32*10^6; % Failure stress; Pa
rho_abs=900; % ABS plastic density; kg/m^3


% Aerofoil Parameters


Ka=0.60; % Coefficient from MIT resource
Ki=0.036; % Coefficient from MIT resource
t_max_ratio=0.10; % Maximum thickness
h_max_ratio=0.051; % Maximum camber
%this is the aero foil data of the airfoil sg6043

% Known Parameters
omega=omega_rpm*2*pi/60; % Angular velocity; rads^-1
TSR=omega*R/Vwind; % Tip Speed Ratio; between 3 and 8
Vt=(1-a1)*Vwind; % Velocity with the induction factor

% Geometry
pan_width=step;

alpha=atan(Vt/omega.*r)-beta;%this is the angle of attack for each segment;
gamma=beta+alpha; % Total angle; in radians


% Force, Moment
for i=1:np
    W(i)=sqrt(Vt^2+(omega*r(i))^2); % Relative wind speed; ms^-1
    L(i)=0.5*rho*W(i)^2*c(i)*Cl(i)*pan_width; % Lift; N
    D(i)=0.5*rho*W(i)^2*c(i)*Cd(i)*pan_width; % Drag; N
    Fthrust(i)=L(i)*cos(gamma(i))+D(i)*sin(gamma(i)); % Aerodynamic force component in the wind direction
    Ftorque(i)=L(i)*sin(gamma(i))-D(i)*cos(gamma(i)); % Aerodynamic force component in the rotatinal direction
end

% In the rotational direction
Msup_to=sum(Ftorque.*r); % Moment at the support
Fsup_to=sum(Ftorque); % Force at the support

Fcut_to=zeros(1,np); % Preallocation
Mcut_to=zeros(1,np); % Preallocation
for i=2:np-1
    Fcut_to(i)=Fsup_to-sum(Ftorque(1:i)); % Shear force at a cut
    cut_dist=r(i)-pan_width/2; % Cut distance
    Mcut_to(i)=Msup_to-Fsup_to*cut_dist-sum(Ftorque(1:i).*r(1:i))+sum(Ftorque(1:i).*cut_dist); %
    M_check(i)=sum(Ftorque(i:np).*(r(i:np)-(r(i-1)+pan_width/2)));
end

% In the wind direction(Thrust)
Msup_th=sum(Fthrust.*r); % Moment at the support
Fsup_th=sum(Fthrust); % Force at the support

Fcut_th=zeros(1,np); % Preallocation
Mcut_th=zeros(1,np); % Preallocation
for i=2:np-1
    Fcut_th(i)=Fsup_th-sum(Fthrust(1:i)); % Shear force at a cut
    cut_dist=r(i)-pan_width/2; % Cut distance
    Mcut_th(i)=Msup_th-Fsup_th*cut_dist-sum(Fthrust(1:i).*r(1:i))+sum(Fthrust(1:i).*cut_dist); %
    Mcheck_th(i)=sum(Fthrust(i:np).*(r(i:np)-(r(i-1)+pan_width/2)));
end

% Stresses
tau=t_max_ratio; % same thing
eps=h_max_ratio; % same thing
A=Ka*c.^2*tau; % cross sectional area of each element of the blade
Itorque=Ki*c.^4*tau*(tau^2+eps^2);% second moment of area

% Stresses due to centrifugal force
for i=1:np-1
Fcentr(i)=rho_abs*omega^2*(pan_width/2)*(A(i)*r(i)+A(np)*r(np)+2*sum((A(i+1:np-1).*r(i+1:np-1))));
end
stress_centr=Fcentr./A(1:np-1);

%stress due to centrifugal force test
if max(stress_centr)>0.8*sigma_max
disp('The centrifugal one doesnt pass');
disp(max(stress_centr));
else
disp('Centrifugal passed');
end

% Stresses due to torque in rotaotional direction test
stress_torque_to=Mcut_to.*(t_max_ratio.*c)./Itorque;
if max(stress_torque_to)>0.8*sigma_max
disp('The torque one in rotational direction doesnt pass');
disp([max(stress_torque_to),'Nm']);
else
disp('Torque in rotational direction passed')
end

%Stresses due to torque in wind direction.
stress_torque_th=Mcut_th.*(t_max_ratio.*c)./Itorque;
if max(stress_torque_th)>0.8*sigma_max
disp('The torque one in wind direction doesnt pass');
else
disp('Torque in wind direction passed');
end

% Stresses due to thrust
fprintf('The maximum limit is  is %2.4f \n',max(0.8*sigma_max));
fprintf('The maximum stress torque in rotational direction is %2.4f Pa\n',max(stress_torque_to));
fprintf('The maximum stress torque in wind direction is %2.4f Pa\n',max(stress_torque_th));
fprintf('The maximum stress due to centerifugal force is %2.4f N\n',max(stress_centr));

function [Cl,Cd,c,Beta,r,step,R]=GetData()

B = 3; %Number of blades
TSR = 5.3; %Theoretical max TSR
a = 1/3; %Induction factor
Vwind = 12; %Wind Speed
R = 0.22; 
m = 20; %Number of sections blade divided into
rho = 1.225;
mu = 1.81 * 10^-5;
r0 = 0.03; %Starting blade radius
cold = zeros([1,101]); %To establish while loop
omega = TSR*Vwind/R; %Angular velocity


%Calculation of radius array
step = (R-r0)/m;

for iiii = 1:m+1
    r(iiii) = r0 + step*(iiii-1);
end

c = (16*pi*R^2)./(9*B*TSR^2*1.03.*r); %Chord variation from given equation
cold = linspace(0,100,m+1); %Initial cold values for comparison in while loop
TSRlocal = omega.*r/Vwind; %Local TSR value at each section
Vrel = sqrt((omega.*r).^2 + (Vwind*(1-a))^2); %Relative wind velocity

while abs(cold-c)>0.0001 %While loop for chord iterationo
    
Re = (rho.*Vrel.*c)./mu; %Calculation of Re at each section

%Interpolation for Re vs Alpha
x = [10000 50000 100000 200000 500000 1000000 2000000];
y = [10 8.75 7 5.5 3.5 2.5 1.5];

polyArray = cubspline(x,y);

%Calculation of optimal AoA at each section using interpolation
for i = 1:m+1
    AoA(i) = cubicEval(x,polyArray,Re(i));
end

AoA = AoA.*pi./180; %Conversion to radians

Gamma = atan((1-a)*R./(TSR.*r)); %Calculation of Gamma which is AoA + Beta (twist)
Beta = Gamma - AoA; 

%Computation of area of each section
for iii = 1:m
    S(iii) = (0.5*(c(iii) + c(iii+1)))*step;
end

%Calculation of average velocity at centre of each section
for iiiii = 1:m
    Vavg(iiiii) = (Vrel(iiiii) + Vrel(iiiii+1))/2;
end

%Reading data file for SG6043 airfoil
fID = fopen('sg6043-2.csv');
filedata = fscanf(fID,'%f %f %f %*f %*f %*f %*f',[3,inf]);
filedata = filedata';

alpha = filedata(:,1);
Clreal = filedata(:,2);
Cdreal = filedata(:,3);

polyArray = cubspline(alpha,Clreal);

AoA = AoA.*180./pi; %Conversion to degrees for interpolation below

%Calculation of optimal Cl at boundary of each section
for iiiiiiii = 1:m+1
    Cl(iiiiiiii) = cubicEval(alpha,polyArray,AoA(iiiiiiii));
end

polyArray = cubspline(alpha,Cdreal);

%Calculation of optimal Cd at boundary of each section
for iiiiiiiii = 1:m+1
    Cd(iiiiiiiii) = cubicEval(alpha,polyArray,AoA(iiiiiiiii));
end

%Calculation of average optimal AoA at centre of each section
for iiiiii = 1:m
    AoAavg(iiiiii) = (AoA(iiiiii) + AoA(iiiiii+1))/2;
end

%Calculation of average optimal Cl at centre of each section
for iiiiiii = 1:m
   Clavg(iiiiiii) = (Cl(iiiiiii) + Cl(iiiiiii+1))/2;
end

cold = c;

F = (2/pi)*acos(exp(-((B*(R-r))./(2.*r.*sin(Gamma))))); %Prandtl Tip Loss Correction

F(21) = F(20); %To ensure chord values result in Re values within range for interpolation

c = (16*pi*R^2.*F)./(9*B*TSR^2.*Cl.*r); %Recalculation of chord values using obtained data (Cl)

end


L = 0.5*rho.*Clavg.*Vavg.^2.*S;
D = 0.5*rho*0.00588*Vavg.^2.*S;
Cp = 0.5926*(1-TSR.*(Cd./Cl));
Ab = R*Vwind/(B*omega);


lift_disp = num2str(2*sum(L));
disp(['Total lift produced is ',lift_disp,' N']);

drag_disp = num2str(2*sum(D));
disp(['Total drag produced is ',drag_disp,' N']);

Area_disp = num2str(2*sum(S));
disp(['Total planform area is ',Area_disp,' m^2'])

Cp_disp = num2str((sum(Cp))/(m+1));
disp(['Average Power Coefficient is ',Cp_disp])

plot(r,c,'k');
grid on

Beta = Beta*180/pi;
Gamma = Gamma*180/pi;
end