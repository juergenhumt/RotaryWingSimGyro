function out = heliDeriv (mu, lamNf, lamI, alfaNf, thtaN, PvL);

global bTl gamRot DvLRes cLsRes aLift thta1 cMI sigma fpA roL vHeli IyHeli mHeliG...
    oM rRot fRotCG hRotCG DvLpGlo cDrgGlo thtaClmb B1LngTrm etaDwn...
    aLiftTail iRigTail cMs mu...
    t21 t22 t23 t24 t25 t26 t31 t32 t33 dt31dmu dt32dmu dt33dmu...
    t11 t12 t13 t14 t15 t16 t17 t18 t19 t11N...
    f2m




hX = 1.0e-11;
coeff716(mu);

dy = eq6(mu, lamNf, thtaN+hX) - eq6(mu, lamNf, thtaN);

dcTsdthta_c_mu = 0.5*aLift*dy/hX; % ./diff(x);
laMu = lamNf/mu;

% cTs from equation 6 NACA716
cTs = 0.5*aLift*eq6(mu, lamNf, thtaN);   
cT  = cTs*sigma;

viBar = lamI/sqrt(0.5*cT);

cX  = (2*mu^2/cT/laMu*(1+laMu^2)^1.5-1);
CXX = (1+ 1/cX)/((1+laMu^2)^0.5);

CNN=(2/aLift+t31*sigma*CXX/2/mu); % eq. (A10) p. 32 NACA2309

% 
cLamMu  = (1 + 1/cX)/sqrt(1+laMu^2);
aN = gamRot*(t11*lamNf+t12*thtaN+t13*thta1)-cMI;
a1 = t14*lamNf+t15*thtaN+t16*thta1;

a2 = mu^2*(t21*lamNf+t22*thtaN+t23*thta1);
b2 = mu^2*(t24*lamNf+t25*thtaN+t26*thta1);

a1Br=2*mu*(4/3*thtaN+lamNf)/(1-0.5*mu^2);
cDx=0.013;
lamP = lamNf+mu*a1;
% hCBram = 0.5*aLift*(0.5*mu*cDx/aLift+a1*thtaN/3+0.75*lamNf*a1-0.5*mu*thtaN*lamNf+0.25*mu*a1^2);
hCBram = 0.25*(mu*cDx*bTl^2+aLift*mu*lamP*bTl*(bTl*lamP+1/3*thtaN*(bTl^2-4.5*mu^2))/(bTl^2+1.5*mu^2));

fprintf('%s%10.5f\n','hC Br   disk            ',hCBram);
fprintf('%s%10.5f\n','4/3*B*thta + lam        ',4/3*bTl*thtaN+lamNf);

alfaD = alfaNf + a1;
fprintf('alfaD    %8.2f°\n',alfaD);
dlamIdmu= (2*mu*thtaN + alfaNf - 4*cTs/aLift/mu)/(1+8*mu/aLift/sigma);

dlamdmu=alfaNf -dlamIdmu;
fprintf('\n%s\n','### a1 Deriv ###'); % 
fprintf('%s%10.5f%s%10.5f%s\n','a1                ',a1,'   a1  Br         ',a1Br,' rad');
% for a1 values an empirical correction proposed by Bramwell
% is added to give better agreement with experiRigMRental values
fprintf('%s%10.5f%s%10.5f%s\n','a1                ',rad2deg(a1)*(1+0.5*mu),'   a1  Br         ',rad2deg(a1Br)*(1+0.5*mu),' rad');

da1dmu=2*(4/3*bTl*thtaN+lamNf)*(bTl^2-1.5*mu^2)/(bTl^2+1.5*mu^2)^2*(1+mu);
fprintf('%s%10.5f%s\n','da1dmu  1            ',da1dmu,' rad');

da1dmu= a1/mu-2*mu/(1-0.5*mu^2)*dlamdmu;   % -- 
fprintf('%s%10.5f%s\n','da1dmu  2            ',da1dmu,' rad');
cAlfa=(bTl^2-0.5*mu^2)*(8*mu+sigma*aLift);

dlamNfdAlfa = 8*mu^2*(bTl+1.5*mu^2)/cAlfa;
da1dalfa   = 16*mu^3/cAlfa;
fprintf('%s%10.5f\n','da1dalfa             ',da1dalfa);
da1dq=-16/(gamRot*(1-0.5*mu^2));
fprintf('%s%10.5f\n','da1dq                ',da1dq);

da1dw=16*mu^2/(1-0.5*mu^2)/(8*mu + aLift*sigma);
fprintf('%s%10.5f\n','da1dw                ',da1dw);


fprintf('\n%s\n','### hc Deriv ###'); % 
dhcdalfa=2/3*bTl*aLift*mu^3*(6*bTl*lamNf+thtaN*(bTl^2-4.5*mu^2))/cAlfa;
fprintf('%s%10.5f\n','dhcdalfa             ',dhcdalfa);
dhcdq = -0.25*aLift*(0.5*lamNf+mu*a1-mu^2*thtaN)*da1dq;
fprintf('%s%10.5f\n','dhcdq                ',dhcdq);
dhcdmu = 0.5*cDx;
fprintf('dhcdmu  %10.5f\n',dhcdmu);
dhcdw = 0.25*aLift/(1+0.25*aLift*lamI/cTs+viBar^4)*(0.5*a1-mu*thtaN+mu*lamP/(1-0.5*mu^2));
fprintf('%s%10.5f\n','dhcdw                ',dhcdw);


nU=0.9;
fprintf('\n%s\n','### b1 Deriv ###'); % 
b1 = gamRot*(t11*lamNf+t18*thtaN+t19*thta1)+t11N*cMI;
b1Br = 4*(mu*aN+ 1.1*nU^0.5*lamNf)/(1-mu^2/2);

fprintf('%s%10.5f%s%10.5f%s\n','b1                ',b1,'   ',rad2deg(b1),'°');
a1S= a1-B1LngTrm;
fprintf('%s%10.5f%s%10.5f%s\n','a1s               ',a1S,'   ',rad2deg(a1S),'°');

% 
fprintf('\n%s\n','### cTs Deriv ###'); % 
dcTsdalfaNf = t31*mu/cos(alfaNf)^2/CNN; % eq. (A11) p. 32 NACA2309
fprintf('%s%10.5f%s\n','dcTsdalfaNf                 ',dcTsdalfaNf,' rad');
dcTsdw = -0.25*aLift/(1+0.25*aLift*(lamI/cTs)+(lamI/sqrt(0.5*cT)));
fprintf('%s%10.5f\n','dcTsdw                   ',dcTsdw);
dcTsdw = -2*aLift*mu/(8*mu+aLift*sigma);
fprintf('%s%10.5f\n','dcTsdw                   ',dcTsdw);

dcTsdthta = t32/CNN; % eq. (A5) p. 31 NACA2309
fprintf('%s%10.5f\n','dcTsdthta  a               ',dcTsdthta);
fprintf('%s%10.5f\n','dcTsdthta  b               ',dcTsdthta_c_mu);

% Coefficients of formula 2 page 13 NACA2309. The coefficients are 
% plotted in figure 1(a) of NACA2309 and are listed on page 34.
k1 = aLift/2*(dt32dmu-t32/t31*dt31dmu) - t32/(mu*CNN);
k2 = dt31dmu/t31+2/(aLift*mu*CNN)+t31*sigma/(mu^2*CNN);
k3 = -(t32/t31)^2*t31/(mu^4*CNN);
k4 =  4*t32/(aLift*t31*mu^4*CNN);
k5 = -4/(aLift^2*t31*mu^4*CNN);

% The lines below give the differences of the values as calculated
% in the example of NACA2309 to those in this program
% d1 = k1*thtaN - (-3.4*0.16);
% fprintf('%s%10.5f\n','d1     ',d1);
% d2 = k2*cTs - (6.2*0.094);
% fprintf('%s%10.5f\n','d2     ',d2);
% d3 = k3*cTs*sigma*thtaN^2 - (-3.1*100*0.094*0.07*0.16^2);
% fprintf('%s%10.5f\n','d3     ',d3);
% d4 = k4*cTs^2*sigma*thtaN - (0.7*1000*0.094^2*0.07*0.16);
% fprintf('%s%10.5f\n','d4     ',d4);
% d5 = k5*cTs^3*sigma - (-3.6*100*0.094^3*0.07);
% fprintf('%s%10.5f\n','d5     ',d5);
% sD = d1+d2+d3+d4+d5;
% fprintf('%s%10.5f\n','Summe d ',sD);
% 

lamP = lamNf+mu*a1;

dtCdMuBr= 2*aLift*mu/(8*mu+aLift*sigma)*(2*mu*thtaN+alfaNf+0.5*cT/mu);
% Formula 2 page 13 NACA2309
dcTsdMu = k1*thtaN + k2*cTs + cT*(k3*thtaN^2 + k4*cTs*thtaN + k5*cTs^2); 
fprintf('%s%10.6f%s%10.6f\n','dcTsdMu                     ',dcTsdMu,'  dtCdMuBr  ',dtCdMuBr);

dcTsdV = dcTsdMu/(oM*rRot);  % eq. (9) p. 16 NACA 2309
fprintf('%s%10.6f\n','dcTsdV   s/m                ',dcTsdV);
fprintf('%s%10.6f\n','dcTsdV   s/ft               ',dcTsdV*f2m);

dcTsdoM = -mu/oM*dcTsdMu; % eq. (10) p. 16 NACA 2309
fprintf('%s%10.6f\n','dcTsdoM                     ',dcTsdoM);

hX = 1.0e-5;
coeff716(mu);
aP = [1 1]; cTsA= [1 1]; tN = [1 1];
[aP(1),cTsA(1),out3]=PvLMuIt(thtaN,mu,PvL);

xh=mu+hX;
coeff716(xh);

pvL(1) = PvL;
options = optimset('MaxFunEvals',5000,'TolFun',1.0e-12,'Display','off');
[pvL(2),fVal] = fsolve(@caPdMu_c_cTsIt,pvL(1)-0.05,options,xh,cTsA(1));
dPvLdMu = diff(pvL)/hX;

[aP(2),cTsA(2),out3]=PvLMuIt(thtaN,xh,pvL(2));
[outS,err]=sprintf('%s','a"');
fprintf('\n%s%s%s\n','### ',outS,' Deriv ###');

daPdMu_c_cTs = (aP(2)-aP(1))/hX;
fprintf('%s%10.6f%s\n','daPdMu   cTs=const          ',daPdMu_c_cTs,'°');
% fprintf('%s%10.6f\n','cTs 1          ',cTs(1));
% fprintf('%s%10.6f\n','cTs 2          ',cTs(2));
% 
% fprintf('%s%10.6f\n','pvL 1          ',pvL(1));
% fprintf('%s%10.6f\n','pvL 2          ',pvL(2));
% 
% fprintf('%s%10.6f%s\n','aP 1          ',aP(1),'°');
% fprintf('%s%10.6f%s\n','aP 2          ',aP(2),s');

% ## '
hX = 1.0e-10;
xh=mu+hX;
coeff716(xh);
[aP(2),cTsA(2),out3]=PvLMuIt(thtaN,xh,PvL);

daPdcTs_c_PvL = diff(aP)./diff(cTsA);
fprintf('%s%10.6f\n','daPdcTs  PvL=const dMu      ',daPdcTs_c_PvL);

PvL = PvL+hX;
coeff716(mu);
[aP(2),cTsA(2),out3]=PvLMuIt(thtaN,mu,PvL);
PvL = PvL-hX;
dPvLdcTs = hX/diff(cTsA);


daPdcTs_c_Mu = diff(aP)./diff(cTsA);
fprintf('%s%10.6f\n','daPdcTs  Mu= const dPvL     ',daPdcTs_c_Mu);

daPdalfaD= daPdcTs_c_Mu*dcTsdalfaNf;
fprintf('%s%10.6f\n','daPdalfaD                   ',deg2rad(daPdalfaD));


fAux  = bTl^3*aLift*thtaN/6/cTs;
daPdq = -16/(gamRot*bTl^4*oM)*0.5*(3-fAux);
fprintf('%s%10.6f\n','fAux                        ',fAux);
fprintf('%s%10.6f\n','3- fAux                     ',(3-fAux));
fprintf('%s%10.6f\n','daPdq                       ',daPdq);


daPdV_c_cTs = deg2rad(daPdMu_c_cTs)/oM/rRot;
fprintf('%s%10.6f\n','daPdV   cTs=const  s/m      ',daPdV_c_cTs);
fprintf('%s%10.6f\n','daPdV   cTs=const  s/ft     ',daPdV_c_cTs*f2m);
xh1=deg2rad(daPdMu_c_cTs); xh2=deg2rad(daPdcTs_c_Mu);
daPdV_c_alfaD = (xh1+xh2*dcTsdMu)/oM/rRot;
fprintf('%s%10.6f\n','daPdV  alfaD=const  r[m]    ',daPdV_c_alfaD);
fprintf('%s%10.6f\n','daPdV  alfaD=const  r[ft]   ',daPdV_c_alfaD*f2m);



% ##
xh=thtaN+hX;
[aP(1),cTs(1),out3]=PvLMuIt(thtaN,mu,pvL(1));
options = optiRigMRset('MaxFunEvals',5000,'TolFun',1.0e-12,'Display','off');
[pvL(2),fVal] = fsolve(@cPvLdThta_c_cTsIt,pvL(1),options,xh,cTs(1));

fprintf('\n%s\n','### PvL Deriv ###');
dPvLdThta_c_cTs = (pvL(2)-pvL(1))/hX;
fprintf('%s%10.6f\n','dPvLdMu                     ',dPvLdMu);
fprintf('%s%10.6f\n','dPvLdcTs                    ',dPvLdcTs);

fprintf('%s%10.6f\n','dPvLdThta   cTs=const       ',dPvLdThta_c_cTs);

fprintf('\n%s\n','### cQs Deriv ###');
dcQsdalfaD = mu*dcTsdalfaNf*(PvL+cTs(1)*dPvLdcTs);
fprintf('%s%10.6f\n','dcQsdalfa                   ',dcQsdalfaD);


fprintf('\n%s\n','### q Deriv ###');
da1dq = -16/gamRot/(1-0.5*mu^2);
fprintf('%s%10.6f\n','da1dq                   ',da1dq);

db1dq = -1/(1+0.5*mu^2);
fprintf('%s%10.6f\n','db1dq                   ',db1dq);

dhcddq = 0.25*aLift*(0.5*lamNf+mu*a1-mu^2*thtaN)*da1dq;
fprintf('%s%10.6f\n','dhcddq                   ',dhcddq);

fprintf('\n%s\n','### Tail Deriv ###');
vBar = vHeli/oM/rRot;

% iRigMR = 0; % Rigging angle of rotor shaft
% iRigTail = 0; % Rigging angle of tail plane 
% etaDwn= 0; % Downwash angle due to main rotor 
fprintf('%s%10.6f\n','etaDwn        ',etaDwn);
epsMH = 2.0*lamNf/vBar; % Downwash angle due to fuselage
fprintf('%s%10.6f\n','epsMH        ',epsMH);

xTail = 1; % Distance from CG to tail plane devided 
           % by rotor Radius
           
jB    = IyHeli/(mHeliG*rRot^2); % Moment of inertia
fprintf('%s%10.6f\n','jB           ',jB);


alfaTL = alfaNf - a1S - iRigMR + iRigTail - (epsMH + etaDwn);
cLT= aLiftTail*(alfaTL - etaDwn);
aT = 0.0044; % Ratio between tail plane and disk area


mUT = -mu*vBar*(cLT + 0.5*aT*(dlamIdmu-lamNf/mu));
fprintf('%s%10.6f\n','mUT          ',mUT);

dlamIdw = 1/(1+0.25*aLift*lamI/cTs+(lamNf/oM/rRot)^2);
mWT = -0.25*mu*vBar*aLiftTail*(1-dlamIdw);
fprintf('%s%10.6f\n','mWT          ',mWT);

mQT = -0.5*mu*aLiftTail*vBar;
fprintf('%s%10.6f\n','mQT          ',mQT);


fprintf('\n%s\n','### Stability Parameter ###');
xU = -(cTs*da1dmu + alfaNf*dcTsdMu + dhcdmu);
fprintf('%s%10.6f\n','xU           ',xU);
xW = -(cTs*da1dw + alfaNf*dcTsdw+dhcdw);
fprintf('%s%10.6f\n','xW           ',xW);
xQ=-(cTs*da1dq + alfaNf*dcTsdw + dhcdq);
fprintf('%s%10.6f\n','xQ           ',xQ);
zU = -dcTsdMu;
fprintf('%s%10.6f\n','zU           ',zU);
zW = -dcTsdalfaNf/mu;
fprintf('%s%10.6f\n','zW           ',zW);


mqf=0; % fuselage pitch contribution
% rotor force leaverage arms
hOne=hRotCG-fRotCG*a1S;
lOne=hRotCG*a1S+fRotCG;

% cMs= moment due to hinge offset
mU = -hOne*xU + lOne*zU + cMs*da1dmu + mqf; % + mUT;
fprintf('%s%10.6f\n','mU           ',mU);

mW = -(fRotCG - hRotCG*a1S)*xW +lOne*zW + cMs*da1dw + mqf; % + mWT;
fprintf('%s%10.6f\n','mW           ',mW);

mQ = -hOne*xQ - 16*cMs/(gamRot*(1-0.5*mu^2)) + mqf; %  + mQT;
fprintf('%s%10.6f\n','mQ           ',mQ);


fprintf('\n%s','The End');
