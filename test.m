clc;
clear;
L=2; Elements=64; Steady=1; Refine=1.5*Steady;
[dx,xspan,dgxspan]=Mesh_Generator(L,Elements,Steady,Refine);
