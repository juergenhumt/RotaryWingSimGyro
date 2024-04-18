clc
clear all
close all


bTl=1.0;
k = 1.5;

% f= bTl^3*aLift/12*thtaN/cTs
f=[0:0.01:3.5];

% xTh = mu*alf/thtaN
xTh = [0.4 0 -0.4];



for m=1:3

for j=1:length(f)

% nnR = (2*(k-1/9) + f*(2/9 + mu*alf/(thtaN*bTl)/3)
  nnR = 2*(k-1/9) + f(j)*(2/9 + xTh(m)/bTl/3);
% dbPdb1Civ= (2*(k-1) + 2*f*(1 + 1.5*mu*alf/(thtaN*bTl)))/nnR;
  y(j)= (2*(k-1) + 2*f(j)*(1 + 1.5*xTh(m)/bTl))/nnR*(0.5/0.29- 0.5*f(j));
end

plot(f,y)
grid
hold on 
end
