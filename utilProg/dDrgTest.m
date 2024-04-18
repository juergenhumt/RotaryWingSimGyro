clc
clear all
close all


dDrgN=     0.0085;     % blade drag coefficient                            :=: dDrgN
dDrg1=    -0.0216;     %  blade drag coefficient                           :=: dDrg1
dDrg2=     0.40;


x=[-0.5:0.01:0.5];

for j=1:length(x)
 f(j)=dDrg2*x(j)^2 + dDrg1*x(j) + dDrgN;
end

plot(x,f,'b')
hold on
grid


dDrgN=     0.0085;     % blade drag coefficient                            :=: dDrgN
dDrg1=    -0.0216;     %  blade drag coefficient                           :=: dDrg1
dDrg2=     0.30;

x=[-0.5:0.01:0.5];

for j=1:length(x)
 f(j)=dDrg2*x(j)^2 + dDrg1*x(j) + dDrgN;
end


plot(x,f,'r')

dDrgN=   0.0105;
dDrg1=  -0.0586;
dDrg2=   0.2700;

x=[-0.5:0.01:0.5];

for j=1:length(x)
 f(j)=dDrg2*x(j)^2 + dDrg1*x(j) + dDrgN;
end


plot(x,f,'g')


dDrgN=   0.0105;
dDrg1=  -0.0286;
dDrg2=   0.2700;

x=[-0.5:0.01:0.5];

for j=1:length(x)
 f(j)=dDrg2*x(j)^2 + dDrg1*x(j) + dDrgN;
end


plot(x,f,'y')
plot(x,f,'-.g')

legend('0.1','0.2','0.3','0.4')


% dDrgN= 0.0085
% ==a=a=== frc func ===a=a==
% xMR    6.500    xcg    6.450    dxCG    0.050    iRigMR    0.0
% cTs 0.07118   cT 0.00427   cLs 6.17909   cL 0.37074
% aN     0.14724°    a1    0.03432°   b1    0.04964°
% Br   aN     8.44°       a1    1.97°      b1    2.84°
% 716  aN     7.54°       a1    1.81°      b1    1.51°
% mmt  aN     6.99°       a1    1.31°      b1    2.84°
% mu  0.148      oM 24.327[1/s]  n 232.3[1/min]   oMR  140.9[m/s]   oMR  462.2[ft/s]
% lamD    0.01668    lamI    0.01436     lamNf    0.01199     mu*a1    0.00470
% A1S   -0.00740°    B1S    0.06844°
% A1S    -0.42°     B1S     3.92°
% alfaNf 10.11°  alfaD 11.92°   alfaS 14.03°   alSm    -2.11°
% thtaN      3.51°       thtaClmb   -15.01      thtaTR     6.78°
% Trim   phi    -0.09°    thta    -0.99°    psi     0.05°
% dDrgN   0.0085  dDrg1  -0.0216  dDrg2   0.4000
% vBx    20.55    vBy    -0.03    vBz     5.13    V    69.48 [fpm]
% vEx    45.88    vEy     0.00    vEz -1080.00    V    47.50 [mph]
% MR    Fx    99.40   Fz -10355.54   Fx/Fz   -0.5°    dxMR    0.0313
% 
% fMR(1)            99.40   fMR(2)            -0.00   fMR(3)        -10355.54
% fTR(1)             0.00   fTR(2)             0.00   fTR(3)             0.00
% fPR(1)           434.20   fPR(2)             0.00   fPR(3)             0.00
% fFus(1)         -692.37   fFus(2)           -3.47   fFus(3)         -971.99
% fFN(1)           -36.79   fFN(2)            20.36   fFN(3)             0.00
% fHS(1)            -0.00   fHS(2)             0.00   fHS(3)            -0.00
% fCG(1)           195.56   fCG(2)           -16.89   fCG(3)         11327.53
% fWN(1)            -0.00   fWN(2)             0.00   fWN(3)             0.00
% fR(1)  -3.1005501e-011   fR(2)   3.0091485e-012   fR(3)  -2.7165675e-011
% 
% Moment Equilibrium
% mMR(1)            -0.00   mMR(2)          -324.45   mMR(3)             0.00
% mTR(1)             0.00   mTR(2)             0.00   mTR(3)             0.00
% mPR(1)             0.00   mPR(2)            -0.00   mPR(3)             0.00
% mFus(1)            0.00   mFus(2)          298.29   mFus(3)           78.49
% mFN(1)             0.00   mFN(2)            26.17   mFN(3)           -78.49
% mHS(1)             0.00   mHS(2)            -0.00   mHS(3)             0.00
% mWN(1)             0.00   mWN(2)            -0.00   mWN(3)             0.00
% mR(1)  -1.8105538e-012   mR(2)   2.2586193e-011   mR(3)  -1.5518253e-011
% ==x=x=== end frc func ===x=x==
% 
% 
% dDrgN=0.0105
% ==a=a=== frc func ===a=a==
% xMR    6.500    xcg    6.450    dxCG    0.050    iRigMR    0.0
% cTs 0.07145   cT 0.00429   cLs 6.19240   cL 0.37154
% aN     0.14858°    a1    0.03534°   b1    0.05003°
% Br   aN     8.51°       a1    2.03°      b1    2.87°
% 716  aN     7.59°       a1    1.87°      b1    1.53°
% mmt  aN     6.99°       a1    1.31°      b1    2.84°
% mu  0.149      oM 24.288[1/s]  n 231.9[1/min]   oMR  140.7[m/s]   oMR  461.5[ft/s]
% lamD    0.01494    lamI    0.01437     lamNf    0.01009     mu*a1    0.00485
% A1S   -0.00810°    B1S    0.07154°
% A1S    -0.46°     B1S     4.10°
% alfaNf  9.35°  alfaD 11.22°   alfaS 13.45°   alSm    -2.23°
% thtaN      3.69°       thtaClmb   -15.01      thtaTR     6.78°
% Trim   phi    -0.09°    thta    -1.57°    psi     0.05°
% dDrgN   0.0105  dDrg1  -0.0516  dDrg2   0.3000
% vBx    20.60    vBy    -0.03    vBz     4.92    V    69.48 [fpm]
% vEx    45.88    vEy     0.00    vEz -1080.00    V    47.50 [mph]
% MR    Fx    91.02   Fz -10360.81   Fx/Fz   -0.5°    dxMR    0.0300
% 
% fMR(1)            91.02   fMR(2)            -0.00   fMR(3)        -10360.81
% fTR(1)             0.00   fTR(2)             0.00   fTR(3)             0.00
% fPR(1)           324.47   fPR(2)             0.00   fPR(3)             0.00
% fFus(1)         -688.40   fFus(2)           -3.44   fFus(3)         -964.17
% fFN(1)           -36.98   fFN(2)            20.45   fFN(3)             0.00
% fHS(1)            -0.00   fHS(2)             0.00   fHS(3)            -0.00
% fCG(1)           309.89   fCG(2)           -17.01   fCG(3)         11324.98
% fWN(1)            -0.00   fWN(2)             0.00   fWN(3)             0.00
% fR(1)  -3.2482796e-011   fR(2)   2.8563818e-012   fR(3)  -2.7883325e-011
% 
% Moment Equilibrium
% mMR(1)            -0.00   mMR(2)          -310.48   mMR(3)             0.00
% mTR(1)             0.00   mTR(2)             0.00   mTR(3)             0.00
% mPR(1)             0.00   mPR(2)            -0.00   mPR(3)             0.00
% mFus(1)            0.00   mFus(2)          284.19   mFus(3)           78.86
% mFN(1)             0.00   mFN(2)            26.30   mFN(3)           -78.86
% mHS(1)             0.00   mHS(2)            -0.00   mHS(3)             0.00
% mWN(1)             0.00   mWN(2)            -0.00   mWN(3)             0.00
% mR(1)  -1.9058461e-012   mR(2)   2.4272336e-011   mR(3)  -1.5177193e-011
% ==x=x=== end frc func ===x=x==

% ==a=a=== frc func ===a=a==
% xMR    6.500    xcg    6.450    dxCG    0.050    iRigMR    0.0
% cTs 0.07162   cT 0.00430   cLs 6.20290   cL 0.37217
% aN     0.14953°    a1    0.03610°   b1    0.05029°
% Br   aN     8.57°       a1    2.07°      b1    2.88°
% 716  aN     7.64°       a1    1.90°      b1    1.54°
% mmt  aN     6.99°       a1    1.32°      b1    2.85°
% mu  0.149      oM 24.263[1/s]  n 231.7[1/min]   oMR  140.5[m/s]   oMR  461.0[ft/s]
% lamD    0.01359    lamI    0.01437     lamNf    0.00863     mu*a1    0.00496
% A1S   -0.00864°    B1S    0.07255°
% A1S    -0.50°     B1S     4.16°
% alfaNf  8.78°  alfaD 10.69°   alfaS 12.94°   alSm    -2.25°
% thtaN      3.82°       thtaClmb   -15.01      thtaTR     6.78°
% Trim   phi    -0.09°    thta    -2.08°    psi     0.05°
% dDrgN   0.0105  dDrg1  -0.0586  dDrg2   0.2700
% vBx    20.64    vBy    -0.03    vBz     4.74    V    69.48 [fpm]
% vEx    45.88    vEy     0.00    vEz -1080.00    V    47.50 [mph]
% MR    Fx    83.52   Fz -10364.56   Fx/Fz   -0.5°    dxMR    0.0287
% 
% fMR(1)            83.52   fMR(2)            -0.00   fMR(3)        -10364.56
% fTR(1)             0.00   fTR(2)             0.00   fTR(3)             0.00
% fPR(1)           226.76   fPR(2)             0.00   fPR(3)             0.00
% fFus(1)         -684.85   fFus(2)           -3.42   fFus(3)         -957.17
% fFN(1)           -37.14   fFN(2)            20.54   fFN(3)             0.00
% fHS(1)            -0.00   fHS(2)             0.00   fHS(3)            -0.00
% fCG(1)           411.71   fCG(2)           -17.12   fCG(3)         11321.74
% fWN(1)            -0.00   fWN(2)             0.00   fWN(3)             0.00
% fR(1)  -3.3686553e-011   fR(2)   2.9487524e-012   fR(3)  -2.8665593e-011
% 
% Moment Equilibrium
% mMR(1)            -0.00   mMR(2)          -297.97   mMR(3)             0.00
% mTR(1)             0.00   mTR(2)             0.00   mTR(3)             0.00
% mPR(1)             0.00   mPR(2)            -0.00   mPR(3)             0.00
% mFus(1)            0.00   mFus(2)          271.56   mFus(3)           79.18
% mFN(1)             0.00   mFN(2)            26.41   mFN(3)           -79.18
% mHS(1)             0.00   mHS(2)            -0.00   mHS(3)             0.00
% mWN(1)             0.00   mWN(2)            -0.00   mWN(3)             0.00
% mR(1)  -1.6199692e-012   mR(2)   2.5557182e-011   mR(3)  -1.4850343e-011
% ==x=x=== end frc func ===x=x==
