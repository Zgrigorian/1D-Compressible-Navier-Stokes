function [Transfer] = Heat_Transfer_Term(states,constants,dx)
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
Transfer=zeros(dgnodes,1);
Transfer(2)=U*Per*dx(1)*(states(Toffset+1)/3+states(Toffset+2)/6-T_surr/2);
for i=2:Elements
    Transfer(2*i-1)=U*Per*dx(i)*(states(Toffset+2*i-1)/3+states(Toffset+2*i)/6-T_surr/2);
    Transfer(2*i)=U*Per*dx(i)*(states(Toffset+2*i-1)/6+states(Toffset+2*i)/3-T_surr/2);
end
Transfer=Transfer(2:end);
end

