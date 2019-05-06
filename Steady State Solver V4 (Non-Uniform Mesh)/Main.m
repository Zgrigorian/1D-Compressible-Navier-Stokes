clc;
%close all;
warning('off','all');
%Mesh/timestepping input variables
L=1; Elements=64; Steady=1; Refine=Steady*1/.7;
%--------------------------------------------------------------------------
%Mesh/timestep initialization
dgnodes=2*Elements;
[dx,xspan,dgxspan]=Mesh_Generator(L,Elements,Steady,Refine);
%--------------------------------------------------------------------------
%Store my Offsets in my vector
uoffset=3*dgnodes; %How many points to offset for velocity part of state
hoffset=0*dgnodes; %How many points to offset for enthalpy part of state
Toffset=2*dgnodes; %How many points to offset for Temperature part of state
rhooffset=1*dgnodes; %How many points to offset for density part of state
Psioffset=4*dgnodes; %How many points to offset for velocity gradient part of state
Gaoffset=5*dgnodes;
Poffset=10*dgnodes;%Pressure is now stored in its own vector
%--------------------------------------------------------------------------
%Relavent Flow/Gas Properties
R_Ideal=8.314; %Ideal gas constant J/mol/K
MW=.029; %Molecular weight (air in this case)
P0=10^5;  %Initial pressure (in Pa)
DP=-1500; %Constant pressure drop (in Pa/m)
R=R_Ideal/MW;  %Ideal gas constant J/g/K
T_surr=290; %Surrounding temperature
Cp=1.0049;  %Heat Capacity (air in this case)
U=10; %Overall Heat transfer Coefficient of pipe wall
Ff=1550; %Constant Frictional Force
Per = 1;  %Perimeter of pipe
A=1;  %Area of pipe
mu=2*10^-6; %Dynamic Viscosity (of air in this case)s
k=2.6*10^-5; %Thermal conductiviity (of air in this case)
u_bc=1; %Velocity Boundary Condition
h_bc=310; %Enthalpy Boundary Condition
%--------------------------------------------------------------------------
poptions=optimset('Display','off');
dogleg='trust-region-dogleg';
trust='trust-region';
leven='levenberg-marquardt';
options=optimset('MaxFunEvals', 500000 ,'MaxIter', 5000,'Display','iter','Algorithm',dogleg);
%--------------------------------------------------------------------------
%Solve for my pressure Profile
P_init=ones(dgnodes,1)*P0;
% P_init=Pressure;
P_func=@(x)Pressure_Corr(x,constants,dx);
[Pressure, Pval]=fsolve(P_func,P_init,poptions);
%--------------------------------------------------------------------------
%Set up my Initial Conditions
% IC=zeros(6*dgnodes,1);
% IC(1:dgnodes)=310;
% IC(2*dgnodes+1:3*dgnodes)=IC(1:dgnodes)/Cp; %Initial Temperature
% IC(dgnodes+1:2*dgnodes)=Pressure./(R*IC(2*dgnodes+1:3*dgnodes)); %Initial Density
% IC(3*dgnodes+1:4*dgnodes)=u_bc;
%IC(4*dgnodes+1:5*dgnodes)=0;
IC=Sol2(end,1:end)';
%--------------------------------------------------------------------------
%Solve for Appropriate Mass Flow Rate
mass_flow=u_bc*IC(dgnodes+1)*A;
constants=[R_Ideal; MW; P0; DP; R; T_surr; Cp; U; Ff; Per; A; mu; k; mass_flow; ...
    uoffset; hoffset; Toffset; rhooffset; Poffset;Elements;dgnodes;Psioffset;Gaoffset]; %Put Everything into a vector
%--------------------------------------------------------------------------
%Use Fsolve to solve the non-linear system of equations
Func=@(x) Solver(x,Pressure,constants,h_bc,dx);
[Sol, fval]=fsolve(Func,IC,options);
%--------------------------------------------------------------------------
%Plot Results
figure(10);
plot(dgxspan,Sol(dgnodes+1:2*dgnodes),'Display', ['Mass Flow rate of ' num2str(mass_flow) ' $\frac{kg}{s}$']);
title('Steady State Density Profile','Interpreter','Latex');
ylabel('Density','Interpreter','Latex');
xlabel('Pipe Coordinate','Interpreter','Latex');
legend('show','Location','best','Interpreter','Latex');
xlim([0,1]);
hold on;
figure(20);
plot(dgxspan,Sol(1:dgnodes),'Display',['Mass Flow rate of ' num2str(mass_flow) ' $\frac{kg}{s}$']);
ylabel('Specific Enthalpy','Interpreter','Latex');
xlabel('Pipe Coordinate','Interpreter','Latex');
tit='Steady State Enthalpy Profile';
title(tit,'Interpreter','Latex');
legend('show','Location','best','Interpreter','Latex');
xlim([0,1]);
hold on;
%legend('Uniform  Mesh','Non Uniform Mesh','Interpreter','Latex');
figure(30);
plot(dgxspan,Sol(uoffset+1:uoffset+dgnodes),'Display',['Mass Flow rate of ' num2str(mass_flow) ' $\frac{kg}{s}$']);
xlabel('Pipe Coordinate','Interpreter','Latex');
ylabel('Velocity','Interpreter','Latex');
tit='Steady State Velocity Profile';
title(tit,'Interpreter','Latex');
legend('show','Location','best','Interpreter','Latex');
xlim([0,1]);
hold on
figure(5);
plot(dgxspan,Sol(hoffset+1:hoffset+dgnodes));
xlabel('Pipe Coordinate','Interpreter','Latex');
ylabel('Enthalpy','Interpreter','Latex');
title('Steady State vs Dynamic Enthalpy Profiles','Interpreter','Latex');
legend('Dynamic Solution at t=5s','Steady State Solution');

figure(6);
plot(dgxspan,Sol(uoffset+1:uoffset+dgnodes));
xlabel('Pipe Coordinate','Interpreter','Latex');
ylabel('Velocity','Interpreter','Latex');
title('Steady State vs Dynamic Velocity Profiles','Interpreter','Latex');
legend('Dynamic Solution at t = 5s','Steady State Solution');

% figure(40);
% plot(dgxspan,Pressure);
% xlabel('Pipe Coordinate','Interpreter','Latex');
% ylabel('Pressure','Interpreter','Latex');
% tit='Steady State Pressure Profile';
% title(tit,'Interpreter','Latex');
% xlim([0,1]);
% %legend('show','Location','best','Interpreter','Latex');
% %legend('Uniform Mesh','Non Uniform Mesh','Interpreter','Latex');
% hold on;
% %Continuity Check
fprintf('The Fsolve Residual is: %2.5s\n',norm(fval));
i=1;
RH=zeros(dgnodes,1);
RH(2*i)=-2*(Sol(rhooffset+2*i)*Sol(uoffset+2*i))/3+...
        Sol(rhooffset+2*i-1)*Sol(uoffset+2*i-1)/3+...
        Sol(rhooffset+2*i-1)*Sol(uoffset+2*i)/6+...
        Sol(rhooffset+2*i)*Sol(uoffset+2*i-1)/6;
for i=2:Elements
    RH(2*i-1)=Sol(rhooffset+2*(i-1))*Sol(uoffset+2*(i-1))-...
             Sol(rhooffset+2*i-1)*Sol(uoffset+2*i-1)/3-...
             Sol(rhooffset+2*i-1)*Sol(uoffset+2*i)/6-...
             Sol(rhooffset+2*i)*Sol(uoffset+2*i-1)/6-...
             Sol(rhooffset+2*i)*Sol(uoffset+2*i)/3;

    RH(2*i)=-2*(Sol(rhooffset+2*i)*Sol(uoffset+2*i))/3+...
        Sol(rhooffset+2*i-1)*Sol(uoffset+2*i-1)/3+...
        Sol(rhooffset+2*i-1)*Sol(uoffset+2*i)/6+...
        Sol(rhooffset+2*i)*Sol(uoffset+2*i-1)/6;
end
fprintf('The Continuity Check Vector has norm: %2.3s\n',norm(RH));

