function [Shear] = Shear_Term(states,constants,dx)
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
Shear=zeros(dgnodes,1);
for i=2:Elements
    Shear(2*i-1)=4/3*A*mu*states(uoffset+2*(i-1))/dx(i-1)*...
                    (states(uoffset+2*(i-1))-states(uoffset+2*(i-1)-1))-...
                    4/3*A*mu*states(uoffset+2*i-1)/dx(i)*...
                    (states(uoffset+2*i)-states(uoffset+2*i-1));
end
Shear=Shear(2:end);
end

