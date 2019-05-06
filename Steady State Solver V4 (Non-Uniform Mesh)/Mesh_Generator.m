function [dx,xspan,dgxspan] = Mesh_Generator(L,Elements,Steady,Refine)
%Mesh/timestep initialization
nodes = Elements + 1; %Number of nodes in regular mesh
dx=zeros(Elements,1); %element length
x_loc=0;
Refine=min(1,Refine);
for i=1:Elements
    if x_loc<=Steady*L
        dx(i)=Refine*L/Elements;
        x_loc=x_loc+dx(i);
    else
        break
    end
end
if i<Elements
    Elements_Left=Elements-(i-1);
    Space_Left=L-x_loc;
    dx(i:Elements)=ones(Elements_Left,1)*Space_Left/Elements_Left;
end
xspan=zeros(nodes,1);
for i=2:nodes
    xspan(i)=xspan(i-1)+dx(i-1);
end
dgxspan = zeros(2*Elements,1); %initializes Dg mesh vector
x_loc=0;
for i = 1:Elements %Generates DG mesh vector
    dgxspan(2*i-1)=x_loc;
    x_loc=x_loc+dx(i);
    dgxspan(2*i) = x_loc;
end
end

