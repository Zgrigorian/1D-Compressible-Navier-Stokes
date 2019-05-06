function [M] = Mass(t,x,constants,Pressure,dx,u_bc,h_bc)
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
M=zeros(4*dgnodes);
M(1:dgnodes-1,:)=Mass_Enth(t,x,constants,dx);
M(dgnodes:2*dgnodes-2,:)=Mass_Cont(t,x,constants,Pressure,dx,u_bc);
end

