function heliDeriv3104(alfaD, thtaN, mu, muZ, oM, cTsReq, lamD, lamI, a1, b1, a1s, vB, anglC, outFile)
%
%
% Copyright 2010/2011 Juergen Humt
% 
% This file is part of RotaryWingSim.
% 
%     RotaryWingSim, is free  software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by the 
%     Free Software Foundation, either version 3 of the License or any later 
%     version.
% 
%     RotaryWingSim is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License along 
%     with RotaryWingSim.  If not, see <http://www.gnu.org/licenses/>.
%
%
%  Version 1.0 of this program was based on report naca-tm-73254 
%
%         --- Version 1.2 ---
%
%
% This routine calculates the stability derivatives of an aircraft using 
% the analytical formulae presented by Bramwell in RaM 3104. Several of
% the input values can be used to compare the values calculated analytically
%
% 
% Input:
%   alfaD   :  disk angle of attack
%   thtaN   :  collective pitch
%   mu      :  advance ratio
%   muZ     :  inflow perpendicular to no feathering plane
%   oM      :  rotor speed
%   cTsReq  :  thrust coefficient
%   lamD    :  disk inflow devided by rotor speed
%   lamI    :  induced inflow devided by rotor speed
%   a1      :  longitudinal flapping coefficient
%   b1      :  lateral flapping coefficint
%   a1s     :  angle between tip path plane and resultant rotor force
%   vB      :  aircraft speed, body coordinates
%   anglC   :  control angles and other values from 
%              main rotor control function
%   outFile :  file pointer
% 

global bTl gamRot aLift dDrgN thta1 cMI sigma aDiskMR fpA vHeli IyHeli mHeliG mHeli...
    iRigMR rRot aBlade xGBlade wBlade eRot wBlade nB zPA zPB xCG hMR...
    dDrgN dDrg1 dDrg2 DvLpGlo cDrgGlo thtaClmb A1LngTrm B1LngTrm etaDwn...
    aLiftHS areaHS lHS hHS iRigHS cMs f2m rhoAir...
    t21 t22 t23 t24 t25 t26 t31 t32 t33 dt31dmu dt32dmu dt33dmu...
    t11 t12 t13 t14 t15 t16 t17 t18 t19 t11N rfrncVal epsGlo

    
fprintf('### heli deriv 3104 ###\n');

coeff716(mu);

oMR = oM*rRot;

cHstrt = 0.25*dDrgN*mu;
fprintf('%s%10.5f\n','cHstrt  ',cHstrt);

doBR = fpA/aBlade;

alfaDF = -0.5*mu^2*doBR*cos(thtaClmb)/cTsReq -cHstrt -thtaClmb;
fprintf('%s%10.3f%s\n','alfaD  ',rad2deg(alfaD),'°');

lamDstrt = -0.5*(fpA/aBlade*mu^3/cTsReq+cTsReq*sigma/mu) -thtaClmb;
fprintf('%s%10.5f\n','lamDstrt  ',lamDstrt);

lamD = muZ-lamI;
fprintf('%s%10.5f%s%10.5f\n','cTs     ',cTsReq,'  cT   ',cTsReq*sigma);

vi    = lamI*oMR;
viBar = lamI/sqrt(0.5*cTsReq*sigma);

ca1Br = (1+0.5*mu);

thtaNBr = thtaN+0.75*thta1;
a1Br2   = 2*mu*(4/3*bTl*thtaNBr+lamD)/(bTl^2+1.5*mu^2);
a1Br    = ca1Br*a1Br2;
a1716 = t14*lamD+t15*thtaN+t16*thta1;
fprintf('%s%10.3f%s%10.3f%s%10.3f%s\n','a1 br D ',rad2deg(a1Br2),'°   a1 Br D mod ',rad2deg(a1Br),'°  a1 716 ',rad2deg(a1716),'°');

b1Br    = 0;
b1716   = gamRot*(t17*lamD+t18*thtaN+t19*thta1);
fprintf('%s%10.3f%s%10.3f%s%10.3f%s\n','b1 br D ',rad2deg(b1Br),'°   b1 Br D mod ',rad2deg(b1Br),'°  b1 716 ',rad2deg(b1716),'°');

lamN=lamD - mu*a1716;
fprintf('%s%10.5f%s%10.5f%s%10.5f%s%10.5f\n','lamD ',lamD,'   lamN ',lamN,'  lamI ',lamI,'  vi ',lamI*oMR/f2m);

xh= bTl^5-0.5*((bTl*mu)^2*(3-5*bTl)+9/4*mu^4);
thtaNloc=1.5*(4/aLift*cTsReq*(bTl^2+1.5*mu^2)-lamN*(bTl^4-0.5*(bTl*mu)^2))/xh;
fprintf('%s%10.3f%s\n','thtaN calc ',rad2deg(thtaNloc),'°');
xh = aLift*lamD*(0.5*a1716-mu*thtaNBr);
hCBr  = 0.25*(mu*dDrgN+xh);
 % hCBr  = 0.25*(mu*dDrgN*bTl^2+ aLift*lamD*mu*bTl*(bTl*lamD+ thtaNBr/3*(bTl^2 - 4.5*mu^2)))/(bTl^2+1.5*mu^2);
hCPr = 0.25*(mu*dDrgN+aLift*mu*lamD*((thtaN*(1/3-1.5*mu^2)+0.5*thta1*(1-1.5*mu^2)+lamD)/(1+1.5*mu^2)));
fprintf('%s%10.5f%s%10.5f\n','hC Pr    ',hCPr,'  hC Br    ',hCBr); 

alfaD = -0.5*fpA/aBlade*mu^2/cTsReq -cHstrt -thtaClmb;
fprintf('%s%10.3f%s\n','alfaD  ',rad2deg(alfaD),'°');

alfaDcalc=atan(lamD/mu+ sigma*cTsReq/(2*bTl^3*mu*(mu^2+lamD^2)^0.5));
fprintf('%s%10.3f%s%10.3f%s%10.3f%s\n','alfaDclc ',rad2deg(alfaDcalc),'°  alfaDinp ',rad2deg(alfaD),'°   alfaD ',rad2deg(alfaD),'°');


qCsBram  = dDrgN*(1+4.7*mu^2)/8-lamN*cTsReq-mu*hCBr;
fprintf('%s%10.5f%s%10.5f\n','qCsBram ',qCsBram,'   qCBram ',qCsBram*sigma);

aPr = atan(hCBr/sqrt(cTsReq^2+hCBr^2));
fprintf('aPr  %10.3f°\n',rad2deg(aPr));

% Start of derivatives from naca-tn-D2655
c1=sigma*aLift/8;
c2=sqrt(mu^2+lamN^2);
cDNTR = (c2/oMR - vi*lamN/c2/oMR^2 + c1/oMR);
alfaNf = alfaD - a1716;


% dvi dthtaN
dvidthtaN=c1*(2/3+mu^2)/cDNTR;
% fprintf('%s%10.5f\n','dvidthtaN ',dvidthtaN);

dvidalfaNf=mu*(lamI^2/c2 + c1*(1-2*thtaNBr*mu*alfaNf))/cDNTR;
% fprintf('%s%10.5f\n','dvidalfaNf ',dvidalfaNf);

dvidV = (c1*(2*thtaNBr*mu*cos(alfaNf)+sin(alfaNf)) - lamI*(mu - lamI*sin(alfaNf))/c2)/(cDNTR*oMR);
% fprintf('%s%10.5f\n','dvidV ',dvidV);

chi = atan(-mu/lamN)+a1716;
c3 = (1-0.5*mu^2);
c4 = (1+0.5*mu^2);

alfaNf = alfaD-a1716;

dhcdmu=0.25*dDrgN*bTl^2;
fprintf('%s%10.5f\n','dhcdmu ',dhcdmu);

dlamdalfaNf=8*mu^2*(bTl^2+1.5*mu^2)/((bTl^2-0.5*mu^2)*(8*mu+sigma*aLift));
fprintf('dlamdalfaNf  %10.5f\n',dlamdalfaNf);

pC=0; zB=1.0e3;  hX=1.0e-8;
cTReq=anglC(15);
alfaD = -tan(anglC(13));  
oM = anglC(8);
xNN2 = sqrt(vB(1)^2 + vB(2)^2);
wNf = alfaD*xNN2;


[cTs1, cT1, lamNf1, lamIj1, kG, jItr] = clcViNf(mu, wNf, oM, thtaN, pC, zB, cTReq);
[cTs2, cT2, lamNf2, lamIj2, kG, jItr] = clcViNf(mu, wNf, oM, thtaN + epsGlo, pC, zB, cTReq);
% xh = cTs2-cTs1;
lamNf = lamNf1;

xh= (lamIj2 -lamIj1)/epsGlo;
fprintf('dlamIthtan %10.5f\n',xh);
xh= (lamNf2 -lamNf1)/epsGlo;
fprintf('dlamdthtan %10.5f\n',xh);


dcTsdalfaD=0.25*bTl^2*aLift*(bTl^2-0.5*mu^2)/(bTl^2+1.5*mu^2)*dlamdalfaNf;

dv1dthtaN = dvidthtaN*tan(0.5*chi) + vi/(1+cos(chi))*(8/3*mu/c3 - mu/oMR*dvidthtaN*(1/c2^2+2/c3));
 fprintf('dv1dthtaN %10.5f   dLamDthtaN\n',dv1dthtaN,dv1dthtaN/oMR);

aN = gamRot/8*(thtaNBr*(1-19/18*mu^2+1.5*mu^4) + 4/3*lamD*c3)/(1+1.5*mu^2);
fprintf('%s%10.3f%s\n','aN   ',rad2deg(aN),'°');

xh= mu^2-mu*lamI*sin(alfaNf)/c2^2 + 2*mu^2/c3 - a1716*tan(alfaNf)*c4/c3;
xh= xh - mu/oMR*dvidalfaNf;
dv1dalfaNf = dvidalfaNf*tan(0.5*chi) + vi/(1+cos(chi))*xh;

xh = lamI*cos(alfaNf)/c2^2 + a1716*cos(alfaNf)/mu*c4/c3;
xh = xh + 2*mu*sin(alfaNf)/c3 - mu*dvidV*(1/c2^2 + 2/c3);
dv1dV = dvidV*tan(0.5*chi) + lamI/(1+cos(chi))*xh;

% clcViNf used to calculate mu and thtaN derivatives 
% since lamN is a function of mu and thtaN
oMR= oM*rRot;

xNN2 = sqrt((vB(1)+hX)^2 + vB(2)^2);
mu2 = xNN2/oMR;
wNf = alfaD*xNN2;
[cTs2, cT2, lamNf2, lamIj, kG, jItr] = clcViNf(mu2, wNf, oM, thtaN, pC, zB, cTReq);

xNN = sqrt(vB(1)^2 + vB(2)^2);
mu1 = xNN/oMR;
wNf = alfaD*xNN;
[cTs1, cT1, lamNf1, lamIj, kG, jItr] = clcViNf(mu1, wNf, oM, thtaN, pC, zB, cTReq);



a1_1=eq2(mu,thtaN,lamNf1);
a1_2=eq2(mu+epsGlo,thtaN,lamNf1);

da1dmu716  = (a1_2 - a1_1)/epsGlo;

fprintf('## a1 deriv ##\n');
xh=4/3*bTl*thtaNBr+lamD;
fprintf('%s%10.5f\n','4/3*B*thta + lam  ',xh);


da1dmu=2*(4/3*bTl*thtaNBr+lamD)*(bTl^2 - 1.5*mu^2)/(bTl^2 +1.5*mu^2)^2*ca1Br;
fprintf('%s%10.5f\n','da1dmu     ',da1dmu);

da1dmu = a1716*c4/c3*cos(alfaNf)/mu + 2*mu/c3*(sin(alfaNf)- dvidV);
fprintf('%s%10.5f\n','da1dmu     ',da1dmu);
fprintf('%s%10.5f\n','da1dmu 716 ',da1dmu716);

da1dthtaN = 2*mu*(4/3*oMR - dvidthtaN)/oMR/c3;
fprintf('%s%10.5f\n','da1dthtaN  ',da1dthtaN);

a1_1=eq2(mu,thtaN,lamNf);
a1_2=eq2(mu,thtaN+epsGlo,lamNf);
da1dthtaN716  = (a1_2 - a1_1)/epsGlo;
fprintf('%s%10.5f\n','da1dthtaN716 ',da1dthtaN716);

a1_1=eq2(mu,thtaN,lamNf);
lamNf2= lamNf+epsGlo;
a1_2=eq2(mu,thtaN,lamNf2);
da1dlam716  = (a1_2 - a1_1)/epsGlo;
fprintf('%s%10.5f\n','da1dlam716 ',da1dlam716);

cNN = (2*mu*(mu^2+lamNf2^2)^0.5);
DvLi = cTReq/cNN;
tanAlfaNf2 = atan(lamNf2/mu+DvLi);
cNN = (2*mu*(mu^2+lamNf^2)^0.5);
DvLi = cTReq/cNN;
tanAlfaNf = atan(lamNf/mu+DvLi);
dAlfa=atan(tanAlfaNf2) - atan(tanAlfaNf);

da1dalfaNf1  = (a1_2 - a1_1)/dAlfa;
fprintf('%s%10.5f\n','da1alfaNf 716 ',da1dalfaNf1);

da1dalfaNf = 2*mu/c3*(mu - dvidalfaNf/oMR) - a1716 * c4/c3*tan(alfaNf);
fprintf('%s%10.5f\n','da1dalfNf ',da1dalfaNf);

da1dw=16*mu^2/(1-0.5*mu^2)/(8*mu + aLift*sigma);
fprintf('%s%10.5f\n','da1dw     ',da1dw);

da1dalfaD = 16*mu^3/(bTl^2-0.5*mu^2)/(8*mu+sigma*aLift)*ca1Br;
fprintf('%s%10.5f\n','da1dalfD  ',da1dalfaD);


da1dq= - 16/gamRot/bTl^4/oM/(1-0.5*mu^2);
fprintf('%s%10.5f\n','da1dq      ',da1dq);

fprintf('## b1 deriv ##\n');
b1Br=4/3*(mu*aN+ 1.1*sqrt(gamRot/16)*lamI)/(1+0.5*mu^2);
fprintf('%s%10.3f%s\n','b1 br D ',rad2deg(b1Br),'°');

db1dthtaN = dv1dthtaN/oMR/c4;
fprintf('%s%10.5f\n','db1dthtaN ',da1dthtaN);

db1dalfaNf =dv1dalfaNf/oMR/c3 - b1Br*mu^2* c4/c3*tan(alfaNf);
fprintf('%s%10.5f\n','db1dalfNf ',da1dalfaNf);

db1dmu = (dv1dV-b1Br*mu*cos(alfaNf))/c3;
fprintf('%s%10.5f\n','db1dmu    ',db1dmu);

dlamIdmu  = (2*mu*thtaN + alfaNf - 4*cTsReq/aLift/mu)/(1+8*mu/aLift/sigma);
dlamNfdmu = 1-dlamIdmu;

db1dmu = gamRot*(lamNf*(0.02091 - mu^2/3) + dlamNfdmu*t17 + (thtaN + 0.75*thta1)*(0.1388 + 0.425*mu^2));
fprintf('%s%10.5f\n','db1dmu    ',db1dmu);

fprintf('## cTs deriv ##\n');
% ----------------
laMu = lamNf/mu;

% cTs from equation 6 NACA716
cTs = cTsReq;   
cT  = cTs*sigma;

cX  = (2*mu^2/cT/laMu*(1+laMu^2)^1.5-1);
CXX = (1+ 1/cX)/((1+laMu^2)^0.5);

CNN=(2/aLift+t31*sigma*CXX/2/mu); % eq. (A10) p. 32 NACA2309


% Coefficients of formula 2 page 13 NACA2309. The coefficients are 
% plotted in figure 1(a) of NACA2309 and are listed on page 34.
k1 = aLift/2*(dt32dmu-t32/t31*dt31dmu) - t32/(mu*CNN);
k2 = dt31dmu/t31+2/(aLift*mu*CNN)+t31*sigma/(mu^2*CNN);
k3 = -(t32/t31)^2*t31/(mu^4*CNN);
k4 =  4*t32/(aLift*t31*mu^4*CNN);
k5 = -4/(aLift^2*t31*mu^4*CNN);

% Formula 2 page 13 NACA2309
dcTsdmu4 = k1*thtaN + k2*cTs + cT*(k3*thtaN^2 + k4*cTs*thtaN + k5*cTs^2); 


% ----------------
c1=2*c1;
dcTsdthtaN = c1*(2/3 + mu^2 - dvidthtaN/oMR);
fprintf('%s%10.5f%s%10.5f\n','dcTsdthtaN ',dcTsdthtaN/sigma,'   dcTdthtaN  ',dcTsdthtaN);

fprintf('%s%10.5f%s%10.5f\n','dcTsdalfaD ',dcTsdalfaD,'   dcTdalfaD ',dcTsdalfaD*sigma);
% fprintf('%s%10.5f%s%10.5f\n','dcTsdalfaD1',dcTsdalfaD1,'   dcTdalfaD1',dcTsdalfaD1*sigma);

dcTsdalfaNf = c1*(mu -2*thtaN*mu^2*tan(alfaNf) - dvidalfaNf/oMR);
fprintf('%s%10.5f%s%10.5f\n','dcTsalfaNf ',dcTsdalfaNf/sigma,'   dcTdalfaNf ',dcTsdalfaNf);

% fprintf('%s%10.5f%s%10.5f\n','dcTsdmu (1)',dcTsdmu1,'   dcTdmu     ',dcTsdmu1*sigma);
dcTsdmu = 2*aLift*mu/(8*mu + aLift*sigma)*(2*thtaN*mu*cos(alfaNf) + sin(alfaNf) + 0.5*cT/mu^2);
fprintf('%s%10.5f%s%10.5f\n','dcTsdmu (2)',dcTsdmu,'   dcTdmu     ',dcTsdmu*sigma);


xh=(1+2.5*(mu-0.17));
dcTsdmu3 = xh*c1*(2*thtaN*mu*cos(alfaNf) + sin(alfaNf) - dvidV)/sigma;
fprintf('%s%10.5f%s%10.5f\n','dcTsdmu (3)',dcTsdmu3,'   dcTdmu     ',dcTsdmu3*sigma);
fprintf('%s%10.5f%s%10.5f\n','dcTsdmu (4)',dcTsdmu4,'   dcTdmu     ',dcTsdmu4*sigma);
dcTsdmu  = dcTsdmu3;

dcTsdw = -0.25*aLift/(1+0.25*aLift*(lamI/cTsReq)+(lamI/sqrt(0.5*(cTsReq*sigma))));
fprintf('%s%10.5f\n','dcTsdw                   ',dcTsdw);
dcTsdw = -2*aLift*mu/(8*mu+aLift*sigma);
fprintf('%s%10.5f\n','dcTsdw                   ',dcTsdw);


fprintf('## cTs deriv autorotation ##\n');
cMM = bTl^2 + 1.5*mu^2;
cYY = sqrt(mu^2 + lamNf^2);
cYY3= cYY^3;
cYY = 2*cYY;
cYY3= 2*cYY3;

cKKK = 1 + 2;


dcTsdalfa = 0.25*bTl^2*(bTl^2 - 0.5*mu^2)/cMM;
dlamdalfa = 0;



fprintf('## cHs deriv ##\n');
dcHsdthtaN = 2/3*a1716 + mu*lamN - dvidthtaN/oMR*(thtaN*mu + 1.5*a1716);
dcHsdthtaN = dcHsdthtaN + da1dthtaN*(1.5*lamN+a1716*mu +2/3*thtaN);
dcHsdthtaN = c1*(dcHsdthtaN + db1dthtaN*mu*lamI/8 + dv1dthtaN*b1Br*mu/8/oMR);
fprintf('%s%10.5f%s%10.5f\n','dcHsdthtaN ',dcHsdthtaN/sigma,'   dcHdthtaN  ',dcHsdthtaN);

dcHsdalfaNf = mu*(thtaN*lamN - dDrgN/aLift - 0.5*a1716^2 -b1Br*lamI/8)*tan(alfaNf) - thtaN*mu^2 +1.5*a1716*mu;
dcHsdalfaNf = dcHsdalfaNf + dvidalfaNf/oMR*(thtaN*mu - 1.5*a1716);
dcHsdalfaNf = c1*(dcHsdalfaNf + da1dalfaNf*(2/3*thtaN + 1.5*lamN + a1716*mu) + db1dalfaNf*mu*lamI/8 + dv1dalfaNf*b1Br*mu/oMR/8);
fprintf('%s%10.5f%s%10.5f\n','dcHsalfaNf ',dcHsdalfaNf/sigma,'   dcHdalfaNf ',dcHsdalfaNf);

dcHsdmu = dDrgN/aLift*cos(alfaNf) - thtaN*(lamN*cos(alfaNf) + mu*sin(alfaNf)) + 1.5*a1716*sin(alfaNf);
dcHsdmu = dcHsdmu + 0.5*a1716^2*cos(alfaNf) + b1Br*lamI*cos(alfaNf)/8;
dcHsdmu = c1*(dcHsdmu + dvidV*(thtaN*mu-1.5*a1716));
fprintf('%s%10.5f%s%10.5f\n','dcHsmu     ',dcHsdmu/sigma,'   dcHdmu     ',dcHsdmu);

dcHsdq = -0.25*aLift*(0.5*lamN+mu*a1716-mu^2*thtaN)*da1dq;
fprintf('%s%10.5f\n','dcHsdq                ',dcHsdq);
dhcdmu = 0.5*dDrgN;
fprintf('%s%10.5f\n','dhcdmu               ',dhcdmu);

viBar = lamI/sqrt(0.5*cTsReq*sigma);
dcHsdw = 0.25*aLift/(1+0.25*aLift*lamI/cTsReq+viBar^4)*(0.5*a1716-mu*thtaN+mu*lamD/(1-0.5*mu^2));
fprintf('%s%10.5f\n','dcHsdw                ',dcHsdw);

dcHsdalfaD=2/3*bTl*aLift*mu^3*(6*bTl*lamD + thtaNBr*(bTl^2- 4.5*mu^2))/(8*mu + sigma*aLift)/(bTl^2-0.5*mu^2);
fprintf('%s%10.5f\n','dcHsdalfaD            ',dcHsdalfaD);


c1=2*c1;
fprintf('## cQs deriv ##\n');
dcQsdthtaN = dv1dthtaN/oMR*0.25*(b1Br- lamI) + 0.25*db1dthtaN*(b1Br+0.5*b1Br*mu^2-lamI) - lamN/3;
dcQsdthtaN = dcQsdthtaN +dvidthtaN/oMR*(thtaN/3 +lamN + 0.5*a1716*mu); 
dcQsdthtaN = c1*(dcQsdthtaN - 0.5*da1dthtaN*(mu*lamN +0.5*a1716 +3/4*a1716*mu^2));
fprintf('%s%10.5f%s%10.5f\n','dcQsdthtaN ',dcQsdthtaN/sigma,'   dcQdthtaN  ',dcHsdthtaN);


dcQsdalfaNf = mu*(lamI - thtaN/3) + mu^2*(0.5*dDrgN/aLift+ 1- (3*a1716^2 -b1Br^2)/8 -0.5*a1716*tan(alfaNf))*tan(alfaNf);
dcQsdalfaNf = dcQsdalfaNf -0.5*a1716*mu^2 - 0.5*a1716*lamI*mu + dvidalfaNf*(thtaN/3+lamN+0.5*a1716*mu)/oMR;
dcQsdalfaNf = dcQsdalfaNf +da1dalfaNf*(0.5*mu*lamI - a1716*(0.25+3/8*mu^2) - 0.5*mu^2*tan(alfaNf));
dcQsdalfaNf = c1*(dcQsdalfaNf + db1dalfaNf*(0.25*lamI-b1Br*(0.25+mu^2/8)) + 0.25*dv1dalfaNf*(b1Br - lamI));
fprintf('%s%10.5f%s%10.5f\n','dcQsdalfaNf ',dcQsdalfaNf/sigma,'   dcQdalfaNf ',dcQsdalfaNf);

dcQsdmu = 0.5*mu*(dDrgN/aLift -0.75 - 0.25*a1716^2)*cos(alfaNf) -(thtaN/3  + lamN +a1716*mu)*sin(alfaNf);
dcQsdmu = dcQsdmu + 0.5*a1716*lamI*cos(alfaNf) + dvidV*(thtaN/3+lamN+0.5*a1716*mu);
dcQsdmu = dcQsdmu - 0.5*da1dmu*(0.5*a1716*(0.5+ 0.75*mu^2) + mu*lamN);
dcQsdmu = c1*(db1dmu + db1dmu*(0.25*lamI-0.25*b1Br*(1 +0.5*mu^2)) + 0.25*dv1dV*(b1Br - lamI));
fprintf('%s%10.5f%s%10.5f\n','dcQsdmu     ',dcQsdmu/sigma,'   dcQdmu ',dcQsdmu);

fprintf('## aPr deriv ##\n');
c5= cos(aPr)^2/(cTsReq*sigma);
daPrdthtaN=(dcHsdthtaN - dcTsdthtaN*tan(aPr))*c5;
fprintf('%s%10.5f\n','daPrdthtaN ',daPrdthtaN);
daPrdalfaNf=(dcHsdalfaNf - dcTsdalfaNf*tan(aPr))*c5;
fprintf('%s%10.5f\n','daPralfaNf ',daPrdalfaNf);
daPrdmu=(dcHsdmu - dcTsdmu*tan(aPr))*c5;
fprintf('%s%10.5f\n','daPrdmu     ',da1dmu);
f1=bTl^3*aLift/6*thtaN/cTsReq;
daPrdq = da1dq*(1.5-0.5*f1);
fprintf('%s%10.5f\n','daPrdq ',daPrdq);

fprintf('\n%s\n','### Tail Deriv ###');

vIx = clcVix(lHS, hHS, vB, anglC);

fprintf('vIx     %10.4f\n',vIx);
epsDwn = vIx/(mu*oMR);

fprintf('%s%10.4f%s%10.4f%s\n','epsDwn   ',epsDwn,'  ',rad2deg(epsDwn),'°');

epsMH = 2.0*lamI/mu; % Downwash angle due to fuselage
fprintf('%s%10.4f%s%10.4f%s\n','epsMH    ',epsMH,'  ',rad2deg(epsMH),'°');

% cMs= moment due to hinge offset
cMs=nB*wBlade*xGBlade*eRot/(2*rhoAir*aBlade*rRot);
fprintf('%s%10.6f\n','cMs    ',cMs);

% Bramwell eq. 6.15 p172
xh=(alfaD+iRigHS-epsDwn);

fprintf('%s%10.4f%s%10.4f%s\n','xh      ',xh,'  ',rad2deg(xh),'°')

vTail = areaHS*lHS;
kBr= 1+0.5*mu^2*vTail*aLiftHS/(cTsReq*hMR + cMs);

cMf=[0 0];

a1sBR= (cMf(1) + hCBr*hMR - cTsReq*xCG-0.5*mu^2*vTail*aLiftHS*(alfaD+iRigHS-epsDwn));
a1sBR=a1s/(kBr*(cTsReq*hMR + cMs));

% Prouty eq. ? p469
% a1s= a1sFnc(thtaN,b1716,a1716,alfaD,mu,cTsReq);    
fprintf('a1sBR     %10.4f  %10.4f°\n',a1sBR,rad2deg(a1sBR));

jB    = IyHeli/(mHeliG*rRot^2); % Moment of inertia
fprintf('%s%10.6f\n','jB           ',jB);

cTs=cTsReq*sigma;

alfaTL = alfaD - a1s - iRigMR + iRigHS - (epsMH + etaDwn);
cLT= aLiftHS*(alfaTL - etaDwn);

if (areaHS > 1.0e-3)
  mUT = -mu*mu*(aLiftHS + 0.5*areaHS*(dlamIdmu-lamI/mu));
  dlamIdw = 1/(1+0.25*aLift*lamI/cTs+lamI/oMR^2);
  mWT = -0.5*mu*mu*aLiftHS*(1-dlamIdw);
  mQT = -0.5*mu*aLiftHS*mu*lHS;
else
  mUT=0;  
  mWT=0;
  mQT=0;
    
end

fprintf('Stability Parameters\n');

fprintf('xU=   -cTs*da1dmu           - alfaD*dcTsdMu      - dhcdmu\n');
fprintf('%s%10.6f%s%10.6f%s%10.6f%s%10.6f%s%10.6f\n','-',cTsReq,'*',da1dmu,' - ',alfaD,'*',dcTsdmu,' - ',dcHsdmu);

xU= -(cTsReq*da1dmu + alfaD*dcTsdmu + dcHsdmu + 2.0*mu*doBR);

fprintf('xW=  -cTs*da1dalfaD        - alfaD*dcTsdalfaD      - dhcdalfaD\n');
fprintf('%s%10.6f%s%10.6f%s%10.6f%s%10.6f%s%10.6f\n','-',cTsReq,'*',da1dalfaD,' - ',alfaD,'*',dcTsdalfaD,' - ',dcHsdalfaD);
xW = -(cTsReq*da1dalfaD + alfaD*dcTsdalfaD + dcHsdalfaD)/mu;

fprintf('%s%10.6f\n','xU     ',xU);  
fprintf('%s%10.6f\n','xW     ',xW);

% fprintf('%s%10.6f\n','xU bra     ',xU*rhoAir*aBlade*oMR^2);
% fprintf('%s%10.6f\n','xU pad     ',xU*mHeli*1000);


xW = -(cTsReq*da1dw + alfaD*dcTsdw + dcHsdw);
fprintf('%s%10.6f\n','xW     ',xW);

fAmr=bTl^3*aLift/6*thtaN/cTsReq;
da1pdq=da1dq*(1.5-0.5*fAmr);
xQ=-(cTsReq*da1pdq*oM  + dcHsdq);
fprintf('%s%10.6f\n','xQ     ',xQ);
zU = -dcTsdmu;
fprintf('%s%10.6f\n','zU     ',zU);
zW = -dcTsdalfaD/mu;
fprintf('%s%10.6f\n','zW     ',zW);

zQ = 16/(gamRot*bTl^4)*dcTsdalfaD;
fprintf('%s%10.6f\n','zQ     ',zQ);

% dlamIdmu= (2*mu*thtaNBr + alfaNf - 4*cTs/aLift/mu)/(1+8*mu/aLift/sigma);
fprintf('dlamIdmu     %10.4f  dlamNfdmu %10.4f°\n',dlamIdmu,dlamNfdmu);



mqf=0; % fuselage pitch contribution
muf=0;
mwf=0;


alfaS = alfaD-a1s;
fprintf('%s%10.3f%s\n','alfaS  ',rad2deg(alfaS),'°');

% rotor force leaver arms
lOne=  xCG       + hMR*alfaS;
hOne= -xCG*alfaS + hMR;

fC= wBlade*xGBlade*oM^2/(rhoAir*aBlade*oMR^2);

mU1 = -hOne*xU + lOne*zU + cMs*da1dmu + muf;
% mU1 = 0.5*fC*eRot*da1dmu - hOne*xU + lOne*zU;
mU = mU1 + mUT;
fprintf('%s%10.6f%s%10.6f%s%10.6f\n','mU     ',mU1,' + mUT ',mUT,'  = ',mU);

mW1 = -hOne*xW +lOne*zW + cMs*da1dw + mwf;
mW  = mW1+ mWT;
fprintf('%s%10.6f%s%10.6f%s%10.6f\n','mW     ',mW1,' + mWT ',mWT,'  = ',mW);

mQ1 = -hOne*xQ + lOne*zQ - 16*cMs/(gamRot*(1-0.5*mu^2)) + mqf;
mQ1 = -hOne*xQ + lOne*zQ +cMs*da1dq + mqf;

mQ = mQ1 +mQT;
fprintf('%s%10.6f%s%10.6f%s%10.6f\n','mQ     ',mQ1,' + mQT ',mQT,'  = ',mQ);


mU=mU/rRot; mW=mW/rRot; mQ=mQ/rRot; 
fprintf('\n   deriv 3104\n');
fprintf('       u                w              q\n');
fprintf('x  %10.6f    %10.6f   %10.6f\n',xU,xW,xQ);
fprintf('z  %10.6f    %10.6f   %10.6f\n',zU,zW,zQ);
fprintf('m  %10.6f    %10.6f   %10.6f\n',mU,mW,mQ);


fprintf(outFile,'\n   deriv 3104\n');
fprintf(outFile,'       u                w              q\n');
fprintf(outFile,'x  %10.6f    %10.6f   %10.6f\n',xU,xW,xQ);
fprintf(outFile,'z  %10.6f    %10.6f   %10.6f\n',zU,zW,zQ);
fprintf(outFile,'m  %10.6f    %10.6f   %10.6f\n',mU,mW,mQ);



fprintf('\n    measured  data\n');
fprintf('       u                w              q\n');
fprintf('x  %10.6f     %10.6f      %10.6f\n',rfrncVal(1),rfrncVal(2),rfrncVal(3));
fprintf('z  %10.6f     %10.6f      %10.6f\n',rfrncVal(4),rfrncVal(5),rfrncVal(6));
fprintf('m  %10.6f     %10.6f      %10.6f\n\n',rfrncVal(7),rfrncVal(8),rfrncVal(9));


% Bramwell 6.14 p167 
B1Br = a1716 + (cMf(1) + hCBr*hMR-cTsReq*xCG)/(cTsReq*hMR+cMs);
fprintf('%s%10.4f%s%10.4f%s\n','B1Br    ',B1Br,'  ',rad2deg(B1Br),'°');

thtaFus=B1Br - a1716 -hCBr/cTsReq -0.5*mu^2*doBR/cTsReq;
fprintf('%s%10.4f%s%10.4f%s\n','thtaFus    ',thtaFus,'  ',rad2deg(thtaFus),'°');

fprintf('### End Deriv 3104 ###\n\n');

fprintf('### Heli Deriv NL ###\n\n');
da1dq= 16/(gamRot*oM*(1-0.5*mu^2));
db1dq= -1/oM/(1+0.5*mu^2);

% dcHdq = 0.5*aLift*sigma*da1dq *(0.25*mu*a1 - 0.5*mu*Ao - 0.25*mu*iC -5/24*lamT)

% xQ_NL= oM*(-cT*da1dq - dcHdq);
fprintf('### End  Deriv NL ###\n\n');


