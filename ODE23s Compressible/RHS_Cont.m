function [RH] = RHS_Cont(t,x,constants,dx)
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
RH=zeros(dgnodes,1);
i=1;
RH(2*i)=-2*(x(rhooffset+2*i)*x(uoffset+2*i))/3+...
        x(rhooffset+2*i-1)*x(uoffset+2*i-1)/3+...
        x(rhooffset+2*i-1)*x(uoffset+2*i)/6+...
        x(rhooffset+2*i)*x(uoffset+2*i-1)/6;
for i=1:Elements
    RH(2*i-1)=x(rhooffset+2*(i-1))*x(uoffset+2*(i-1))-...
             x(rhooffset+2*i-1)*x(uoffset+2*i-1)/3-...
             x(rhooffset+2*i-1)*x(uoffset+2*i)/6-...
             x(rhooffset+2*i)*x(uoffset+2*i-1)/6-...
             x(rhooffset+2*i)*x(uoffset+2*i)/3;

    RH(2*i)=-2*(x(rhooffset+2*i)*x(uoffset+2*i))/3+...
        x(rhooffset+2*i-1)*x(uoffset+2*i-1)/3+...
        x(rhooffset+2*i-1)*x(uoffset+2*i)/6+...
        x(rhooffset+2*i)*x(uoffset+2*i-1)/6;
end
RH=RH(2:end);
end

