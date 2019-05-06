function [F] = Solver_Weak(states,Pressure,constants,h_bc,dx)
%Builds my RHS for my 1st continuity Equation
R_Ideal=constants(1);
MW=constants(2);
P0=constants(3);
DP=constants(4);
R=constants(5);
T_surr=constants(6);
Cp = constants(7);
U =constants(8);
Ff=constants(9);
Per = constants(10);
A = constants(11);
mu = constants(12);
k = constants(13);
m_flow=constants(14);
uoffset= constants(15);
hoffset= constants(16);
Toffset = constants(17);
rhooffset = constants(18);
Poffset = constants(19);
Elements=constants(20);
dgnodes=constants(21);
F=zeros(4*dgnodes,1);
F(1:dgnodes)=Continuity_EQ(states,constants); %Continuity Equations
F(dgnodes+1:2*dgnodes)=Density_Corr(states,Pressure,constants); %Density Correlations
F(2*dgnodes+1:3*dgnodes)=Temp_Corr(states,constants); %Temperature Correlations
F(3*dgnodes+1:4*dgnodes)=Enthalpy_Term(states,constants,h_bc)+Pressure_Term(states,constants,dx)-...
     Friction_Term(states,constants,dx)+Heat_Transfer_Term(states,constants,dx)+Viscous_Term(states,constants,dx);...
     +Heat_Conductance_Term(states,constants,dx)-Shear_Term(states,constants,dx);
% F(4*dgnodes+1:5*dgnodes)=Psi_EQ(states,constants,dx);
% F(5*dgnodes+1:6*dgnodes)=Gamma_EQ(states,constants,dx);
end

