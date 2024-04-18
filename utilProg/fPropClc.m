clc
clear all
close all

vRange=[1:1:200];
rProp= 6*0.3048;
FOM=0.78;

options = optimset('MaxFunEvals',5000,'TolFun',1.0e-10,'Diagnostics','off');

xS=0;
n=length(vRange);
pAvail=75000;

for j=1:length(vRange)
 vR=vRange(j);
 etaR(j)=fsolve(@fProp,xS,options,vR,rProp,pAvail);
 thrst(j)=etaR(j)*pAvail/vR;
 vI(j)= vR*(FOM/etaR(j)-1.0);
end


rProp=4*0.3048;

for j=1:length(vRange)
 vR=vRange(j);
 etaR(j)=fsolve(@fProp,xS,options,vR,rProp,pAvail);
 thrst2(j)=etaR(j)*pAvail/vR;
  vI2(j)= vR*(FOM/etaR(j)-1.0);
end


plot(vRange,thrst);
hold on
title('Thrust vs Speed');
grid
plot(vRange,thrst2,'g');
legend('R=1.2','R=0.8');

figure(2)
plot(vRange,vI);
hold on
title('Induced Velocity vs Speed');
grid
plot(vRange,vI2,'g');
legend('R= 6''','R= 4''');
