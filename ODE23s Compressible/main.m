clc;
close all;
L=1; Elements=64; Steady=1; Refine=4*Steady;
t0=0; tf=5;
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
Poffset=10*dgnodes;%Pressure is now stored in its own vector
Psioffset=4*dgnodes; %Velocity Gradient offset vector
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
mass_flow=1; %Mass flow rate (used in steady state solver)
h_bc=270; %Enthalpy Boundary Condition (Dirchlet)
u_bc=1; %Velocity Boundary Condition (Dirchlet)
constants=[R_Ideal; MW; P0; DP; R; T_surr; Cp; U; Ff; Per; A; mu; k; mass_flow; ...
    uoffset; hoffset; Toffset; rhooffset; Poffset;Elements;dgnodes;Psioffset]; %Put Everything into a vector
options=optimset('MaxFunEvals', 50000 ,'MaxIter', 5000,'Display','off'); %Define max Iterations
%--------------------------------------------------------------------------
%Solve for my pressure Profile
P_init=ones(dgnodes,1)*P0; %Initial Guess for P
P_func=@(x)Pressure_Corr(x,constants,dx,0); %Define Auntonomous Pressure Function
[Pressure, Pval]=fsolve(P_func,P_init,options); %Use Fsolve to solve for Pressure
%--------------------------------------------------------------------------
%Set up my Initial Conditions
IC=zeros(4*dgnodes,1);
IC(1:dgnodes)=h_bc; %Initial Enthalpy
% IC(2:dgnodes)=300;
IC(2*dgnodes+1:3*dgnodes)=IC(1:dgnodes)/Cp; %Initial Temperature
IC(dgnodes+1:2*dgnodes)=Pressure./(R*IC(2*dgnodes+1:3*dgnodes)); %Initial Density
IC(3*dgnodes+1:4*dgnodes)=u_bc; %Initial Velocity
%--------------------------------------------------------------------------
%Initialize my Mass Matrix
M=@(t,x) Mass(t,x,constants,Pressure,dx,u_bc,h_bc);
options=odeset('Mass',M,'MassSingular','yes');
Func=@(t,x) RHS(t,x,constants,Pressure,dx,h_bc,u_bc);
tspan=[t0 tf];
[t,Sol]=ode23t(Func,tspan,IC,options);
test= ode23t(Func,tspan,IC,options);
%--------------------------------------------------------------------------
% figure(1);
% surf(dgxspan,t,Sol(:,1:dgnodes));
% xlabel('Pipe Position','interpreter','Latex');
% ylabel('Time','interpreter','Latex');
% zlabel('Enthalpy','interpreter','Latex');
% title('Enthalpy Profile','interpreter','Latex');
% figure(2);
% surf(dgxspan,t,Sol(:,dgnodes+1:2*dgnodes));
% xlabel('Pipe Position','interpreter','Latex');
% ylabel('Time','interpreter','Latex');
% zlabel('Density','interpreter','Latex');
% title('Density Profile','interpreter','Latex');
% figure(3);
% surf(dgxspan,t,Sol(:,2*dgnodes+1:3*dgnodes));
% xlabel('Pipe Position','interpreter','Latex');
% ylabel('Time','interpreter','Latex');
% zlabel('Temperature','interpreter','Latex');
% title('Temperature Profile','interpreter','Latex');
% figure(4);
% surf(dgxspan,t,Sol(:,3*dgnodes+1:4*dgnodes));
% xlabel('Pipe Position','interpreter','Latex');
% ylabel('Time','interpreter','Latex');
% zlabel('Velocity','interpreter','Latex');
% title('Velocity Profile','interpreter','Latex');
figure(500);
plot(dgxspan,Sol(end,1:dgnodes));
title(['Final Enthalpy at time=' num2str(t(end)) ' with ' num2str(Elements) ' Elements'],'interpreter','Latex');
xlabel('Pipe position','interpreter','Latex');
ylabel('Enthalpy','interpreter','Latex');
hold on;
figure(600);
plot(dgxspan,Sol(end,3*dgnodes+1:4*dgnodes))%,'DisplayName',['t=',num2str(tf)]);
title(['Final Velocity Profile at time=' num2str(t(end)) ' with ' num2str(Elements) ' Elements'],'interpreter','Latex' );
% title(['Final Velocity Profile with ' num2str(Elements) ' Elements'],'interpreter','latex');
xlabel('Pipe Position','interpreter','Latex');
ylabel('Velocity','interpreter','Latex');
hold on;
% % figure(7);
% plot(dgxspan,Sol(1,3*dgnodes+1:4*dgnodes))
% title('Initial Velocity Profile','interpreter','Latex');
% xlabel('Pipe Position','interpreter','Latex');
% ylabel('Velocity','interpreter','Latex');
% getfield(test.stats,'nfevals')
Sol2=Sol;
test.stats.nfevals
disp(t(end));
close all;