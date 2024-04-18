clc

gamRot= 6;
aLift = 5.7;
sigma = 0.06;
cT = 0.00055;
oMR = 157;

R=1;
Z=0.5;

c1 = 1 - (R/(4*Z))^2;

lamI = sqrt(0.5*cT);
xh = 4/aLift*cT/sigma;

thtaN = 4/aLift*cT/sigma + lamI;
thtaNdeg = rad2deg(thtaN)

aNBr  = 0.125*gamRot*( thtaN + 4/3*lamI)
aNBrDeg = rad2deg(aNBr)

lamI = 0.75*lamI;

thtaN = 4/aLift*cT/sigma + lamI;
thtaNdeg = rad2deg(thtaN)

aNBr  = 0.125*gamRot*( thtaN + 4/3*lamI)
aNBrDegIG = rad2deg(aNBr)

xh = 1-aNBrDeg/aNBrDegIG


aN=2/3*12.1*((0.0055/0.06)/5.7 - 1.5*9.81*rRot/oMR^2)
aNdeg = rad2deg(aN)