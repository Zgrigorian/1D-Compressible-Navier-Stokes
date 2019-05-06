function [F] = Psi_EQ(states,constants,dx)
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
Psioffset=constants(22);
F=zeros(dgnodes,1);
for i=1:Elements-1
    F(2*i-1)=states(Psioffset+2*i-1)+(states(uoffset+2*i-1)-3*states(uoffset+2*i)+1.95*states(uoffset+2*i+1))/(dx(i));
    F(2*i)=states(Psioffset+2*i)+(states(uoffset+2*i-1)+3*states(uoffset+2*i)-states(uoffset+2*i+1))/(dx(i));
end
i=Elements;
    F(2*i-1)=states(Psioffset+2*i-1)+10;
    F(2*i)=states(Psioffset+2*i)+10;
% for i=1:Elements-1
%     F(2*i-1)=states(Psioffset+2*i-1)+(states(uoffset+2*i-1)-3*states(uoffset+2*i)+2*states(uoffset+2*i+1))/dx(i);
%     F(2*i)=states(Psioffset+2*i)+(states(uoffset+2*i-1)+3*states(uoffset+2*i)-4*states(uoffset+2*i+1))/dx(i);
% end
 i=Elements;
     F(2*i-1)=states(Psioffset+2*i-1)+(states(uoffset+2*i-1)-states(uoffset+2*i))/dx(i);
     F(2*i)=states(Psioffset+2*i)+(states(uoffset+2*i-1)-states(uoffset+2*i))/dx(i);
end