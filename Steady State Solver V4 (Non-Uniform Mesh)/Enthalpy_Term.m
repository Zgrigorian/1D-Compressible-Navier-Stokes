function [Enth] = Enthalpy_Term(states,constants,h_bc)
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
Enth=zeros(dgnodes,1);
Enth(1)=states(hoffset+1)-h_bc;
Enth(2)=m_flow*A*states(hoffset+2)-m_flow*A/2*(states(hoffset+1)+states(hoffset+2));
for i=2:Elements
    Enth(2*i-1)=-m_flow*A*states(hoffset+2*(i-1))+...
                 m_flow*A/2*(states(hoffset+2*i)+states(hoffset+2*i-1));
    Enth(2*i)=m_flow*A*states(hoffset+2*i)-...
                 m_flow*A/2*(states(hoffset+2*i)+states(hoffset+2*i-1));
end

