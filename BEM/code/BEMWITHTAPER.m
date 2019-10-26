%This code is to use blade element theory, together with optimum Cl, Cd and
%alpha values for an airfoil to optimise for a, find optimal chord length
%and twist angle.

%housekeeping
clear
clc

%Ask for inputs
%Design characteristics and testing parameters
B = 3;%input('Number of blades:');
R = 0.22;%input('Blade Radius:');
r0 = 0.036;%input('Blade starting radius:');
vwind = 12;%input('Wind speed:');
TSR = 4*pi/B;%input('TSR:');
omega = TSR*vwind/R; %Angular velocity

%Airfoil characteristics
alpha = 8.5;%input('Optimum alpha');
Cl = 1.431;%input('Optimum Cl');
Cd = 0.036;%input('Optimum Cd');
rho = 1.225; %density of air

%initialize blade sections (31)
dr=(r0:((R)-r0)/20:R); 

%start BEM algo to find a, chord and Cl over the span of the blade
for i =1:3
    difference1  = 2; %initialize difference
    %Intialize a and aprime
    a = 0;
    aprime = 0;
    while difference1 > 0.01
        r = dr(i);
        vtan = r * omega; %tangential velocity at that blade section
        phi = atan(((1-a)/(1+aprime))*(vwind/vtan)); %Calculating incident angle
        
        %calculating Cx and Cy
        Cx = Cl*sin(phi) - Cd*cos(phi);
        Cy = Cl*cos(phi) + Cd*sin(phi);
        
        %calculating tapered chord length c and solidity
        if i == 1
            c = 0.035;
        elseif i == 2
            c = 0.038;
        elseif i == 3
            c = 0.040;
        end
        sigma = c*B/(2*pi*r); %solidity
        
        %calculate tip loss factor F
        F = (2/pi)*acos(exp((-B/2)*((R-r)/(r*sin(phi)))));
        ac = 0.2; %initializing ac
        aprimedummy = ((4*F*sin(phi)*cos(phi))/(sigma*Cx))-1;
        aprimenew = 1/aprimedummy;
        if a <= 0.2
           anew = 1/(((4*F*sin(phi)^2)/(sigma*Cy))+1);
        elseif a > 0.2
            K = (4*F*sin(phi)^2)/(sigma*Cy);
            anew = 0.5*(2+K*(1-2*ac)-sqrt((K*(1-2*ac)+2)^(2)+4*(K*ac*ac-1)));
        end
        
        %calculating difference to solve for a
            difference1 = abs(a-anew);
            a = anew;
            aprime = aprimenew;
    end
    %inserting values into array
    dc(i) = c; %inserting chord length into array
    beta(i) = (phi*180/pi)-alpha; %inserting twist angle into array
    dCx(i) = Cx;
    dCy(i) = Cy;
    w(i) = sqrt(((1-a)*vwind)^2 + ((1+aprime)*vtan)^2);
    dUcalc(i) =  (vwind*(1-a)*omega*r*(1+aprime))/(sin(phi)*cos(phi));%value to multiply to calculate dU
    dTcalc(i) = ((vwind*(1-a))^2)/(sin(phi)^2); %value to multiply for each section to calculate dT
end
for i = 4:length(dr)
    difference1  = 2; %initialize difference
    %Intialize a and aprime
    a = 0;
    aprime = 0;
    while difference1 > 0.01
        r = dr(i);
        vtan = r * omega; %tangential velocity at that blade section
        phi = atan(((1-a)/(1+aprime))*(vwind/vtan)); %Calculating incident angle
        
        %calculating Cx and Cy
        Cx = Cl*sin(phi) - Cd*cos(phi);
        Cy = Cl*cos(phi) + Cd*sin(phi);
        
        %calculating chord length c and solidity
        c = (1/B)*(16*pi*r/Cl)*(sin((1/3)*atan(R/(TSR*r))))^2;
        sigma = c*B/(2*pi*r); %solidity
        
        %calculate tip loss factor F
        F = (2/pi)*acos(exp((-B/2)*((R-r)/(r*sin(phi)))));
        ac = 0.2; %initializing ac
        aprimedummy = ((4*F*sin(phi)*cos(phi))/(sigma*Cx))-1;
        aprimenew = 1/aprimedummy;
        if a <= 0.2
           anew = 1/(((4*F*sin(phi)^2)/(sigma*Cy))+1);
        elseif a > 0.2
            K = (4*F*sin(phi)^2)/(sigma*Cy);
            anew = 0.5*(2+K*(1-2*ac)-sqrt((K*(1-2*ac)+2)^(2)+4*(K*ac*ac-1)));
        end
        
        %calculating difference to solve for a
            difference1 = abs(a-anew);
            a = anew;
            aprime = aprimenew;
    end
    %inserting values into array
    dc(i) = c; %inserting chord length into array
    beta(i) = (phi*180/pi)-alpha; %inserting twist angle into array
    dCx(i) = Cx;
    dCy(i) = Cy;
    w(i) = sqrt(((1-a)*vwind)^2 + ((1+aprime)*vtan)^2);
    dUcalc(i) =  (vwind*(1-a)*omega*r*(1+aprime))/(sin(phi)*cos(phi));%value to multiply to calculate dU
    dTcalc(i) = ((vwind*(1-a))^2)/(sin(phi)^2); %value to multiply for each section to calculate dT

end

%calculating forces
dU = 0.5*rho.*dUcalc.*dc.*dCx; %Tangential (Vertical force)
dM = 0.5*rho.*dUcalc.*dc.*dCx.*dr; %Torque
dT = 0.5*rho.*dTcalc.*dc.*dCy; %Axial Force
plot(dr,dc,'r-');

Utotal = trapz(dr(1:20),dU(1:20));
Mtotal = trapz(dr(1:20),dM(1:20));
Ttotal = trapz(dr(1:20),dT(1:20));
Ptotal = (omega*B)*trapz(dr(1:20),dU(1:20));
Cp = (0.5*rho*(pi*R^2)*vwind^3)/Ptotal;