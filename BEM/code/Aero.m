
clc
clear


B = 3; %Number of blades
TSR = 4*pi/B; %Theoretical max TSR
a = 1/3; %Induction factor
Vwind = 12; %Wind Speed
R=0.22;
m = input('input the number of the slices of the blades'); %Number of sections blade divided into
rho = 1.225;
mu = 1.81 * 10^-5;
r0 = 0.036; %Starting blade radius
cold = zeros([1,m+1]); %To establish while loop
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

F(m+1) = F(m); %To ensure chord values result in Re values within range for interpolation

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
