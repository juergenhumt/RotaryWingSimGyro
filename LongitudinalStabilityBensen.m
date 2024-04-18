
clc
clear
close All


global dDrgN dDrg1 dDrg2 thta1 bTl gamRot aLift mHeliG...
   sigma rhoAir iRigMR...
   t11 t12 t13 t14 t15 t16 t17 t18 t19 t11N...
   t21 t22 t23 t24 t25 t26 t31 t32 t33 dt31dmu dt32dmu dt33dmu...
   t41 t42 t43 t44 t45 t46 t47 t48 t49 t41N...
   t51 t52 t53 t54 t55 t56 t57 t58 t59 t51N t511 t512 t513 t514...
   t61 t62 t63 t64 t65 t66 t67 t68 t69 t61N t611 t612 t613 t614...    
   t71 t72 t73;

f2m=0.3048;
rhoAir=1.225;
rRot=3.1;
iRigMR=0;
sigma = 0.07;
lbForce= 0.435*9.81;

xCG  =0.1;
hMR  =0.25;

thtaN=deg2rad(2.5);
dDrgN= 0.0165;
dDrg1= -0.21;
dDrg2= 0.5;

nRot=406;
oM=pi*nRot/30;

mHeliG=330*9.81;
cT_const=mHeliG/(rhoAir*pi*rRot*rRot*(oM*rRot)^2);

lamHovConst=sqrt(0.5*cT_const);

aLift=5.7;
thta1=0;

bTl=0.97;
gamRot=10;

muAr = [0.1:0.005:0.2];
lamAr=muAr;
oM   =muAr;
V    =muAr;
doMdV =muAr;
doMdV2=muAr;


alNf =muAr;
xhAr=muAr;


dcTdoM = muAr;
dcTdalf= muAr;

dTdV=muAr;
dTdV_Qn=muAr;


dTdAlf=muAr;


cDj=0;

rRotBs=rRot;

eps1N=10*eps;

for j=1:length(muAr)
    
  mu=muAr(j);
  muPe=mu + eps1N;
  [cTs, cT, lamNf, lamI, alfaD, alfaDeg, alfaNf]=clcRotState2(thtaN,mu);
  [cTsPe, cTPe, lamNfPe, lamIPe, alfaDPe, alfaDegPe, alfaNfPe]=clcRotState2(thtaN,muPe);
  
  
  
  cTsAr(j)=cTs;
  alfaDegAr(j)=alfaDeg;
  
  a1=eq2(mu,lamNf,thtaN);
  b1=eq3(mu,lamNf,thtaN);
  
  % lamNf=lamNf*0.03/0.0262;
  
  lamAr(j)= lamNf;
  lam_Hov(j)= lamNf/lamHovConst;
  mu_lamHov(j)=mu/lamHovConst;
  
  oMR = sqrt(mHeliG/(cT*pi*rRotBs^2*rhoAir));
  
  
  K_1 = rhoAir*pi*rRotBs^2*oMR^2;

  dTdVPe(j) = K_1*(cTPe - cT)/eps1N*cos(alfaNf)*oMR;
  dTdAlfPe(j) = K_1*(cTPe - cT)/(alfaNfPe-alfaNf);

  
  
  oM(j)= oMR/rRotBs;
  vi = lamI*oMR;
  
  oMRPe = sqrt(mHeliG/(cTPe*pi*rRotBs^2*rhoAir));
  viPe=lamIPe*oMRPe;
  
  xh=viPe-vi;
  xh=oMRPe - oMR;
  
  V(j)=muAr(j)/cos(alfaNf)*oM(j)*rRotBs;
  
  
  alNf(j)=alfaNf;

  vHeli=mu*oMR/cos(alfaNf);

  cH(j)   = 0.5*sigma*aLift*cH2sa(mu,lamNf,thtaN);
  res=heliDeriv2655(mu, alfaNf, lamNf, lamI, thtaN, a1, cT, cH(j), oMR, vHeli);
  
  
% for k=1:length(res)
%     fprintf('%10.5f\n',res(k));
% end
% fprintf('dHdalf  %10.5f    dcHdalf %10.5f\n',K_1*res(15),res(15));

dHdalf_26(j)=res(15);

%   [ aN a1 b1 chi real(da1dthtaN) real(da1dalfaNf) real(da1dmu)....
%     real(db1dthtaN) real(db1dalfaNf) real(db1dmu)...
%     real(dcTdthtaN) real(dcTdalfaNf) real(dcTdmu)...
%     dcHdthtaN dcHdalfaNf dcHdmu...
%     dcQdthtaN dcQdalfaNf dcQdmu...
%     daPrdthtaN daPrdalfaNf daPrdmu dvidV]=res;
%   
  t41m=(5/4 + gamRot^2*bTl^8/1296);
  t41x = bTl^2/2 + t41m*mu^2;

  t42m= (8*bTl/3 + gamRot^2*bTl^9);
  t42x = bTl^3/3 + t42m*mu^2;

  t52x = 1/3;
  
  t55m = (-1/4+1/bTl^2+1/2/bTl^4+gamRot^2*(1/162*bTl^4-1/81*bTl^5+1/144*bTl^6));
  t55x = 0.5 +  t55m*mu^2;

  t56m = (4/3/bTl+4/3/bTl^3+gamRot^2*(1/108*bTl^5-1/54*bTl^6+1/96*bTl^7));
  t56x = 2/3+ t56m*mu^2;


  t58m = (1/4+8/9/bTl^2+gamRot^2*(1/288*bTl^6-1/144*bTl^7+1/256*bTl^8));

  dlamdV=sin(alfaNf)/oM(j)/rRotBs;
  
  mu2V = 2*mu^2/V(j);

  xNM = aLift*(t41m*mu2V*lamNf^2 + t41x*2*lamNf*dlamdV +...
           (t42m*mu2V*lamNf + t42x*dlamdV)*thtaN +...
           (8*bTl^2/9 + gamRot^2*bTl^10/1080)*mu2V*thtaN^2);
       
       
  xNM = xNM - (dDrgN*0.25*mu2V +...
     dDrg1*(1/3*dlamdV + 0.25*mu2V*thtaN) +...
     dDrg2*((t55m*mu2V*lamNf^2 + t55x*2*lamNf*dlamdV) +...
             (t56m*mu2V*lamNf +  t56x*dlamdV)*thtaN +...
            t58m*mu2V*thtaN^2));


  dlamdoM= -lamNf/oM(j);
  mu2oM  = -mu^2/oM(j);
 

  xDN = aLift*(t41m*mu2oM*lamNf^2 + t41x*2*lamNf*dlamdoM +...
           (t42m*mu2oM*lamNf + t42x*dlamdoM)*thtaN +...
           (8*bTl^2/9 + gamRot^2*bTl^10/1080)*mu2oM*thtaN^2);
       
       
  xDN = xDN - (dDrgN*0.25*mu2oM +...
     dDrg1*(1/3*dlamdoM + 0.25*mu2oM*thtaN) +...
     dDrg2*((t55m*mu2oM*lamNf^2 + t55x*2*lamNf*dlamdoM) +...
             (t56m*mu2oM*lamNf +  t56x*dlamdoM)*thtaN +...
            t58m*mu2oM*thtaN^2));

  doMdV(j)=-xNM/xDN;

  eps=1.0e-10;
  
  lamNf1=((V(j))*sin(alfaNf)-vi)/oM(j)/rRotBs;
  mu1 = (V(j))*cos(alfaNf)/oM(j)/rRotBs;
  
  lamNf2= ((V(j)+eps)*sin(alfaNf)-vi)/oM(j)/rRotBs;
  mu2   = (V(j)+eps)*cos(alfaNf)/oM(j)/rRotBs;
  
  
%   dcTdV(j)= 0.5*aLift*sigma*(mu2V*( eq6(mu2, lamNf1, thtaN) -  eq6(mu1, lamNf1, thtaN))/eps +...
%       dlamdV*(eq6(mu1, lamNf2, thtaN) -  eq6(mu1, lamNf1, thtaN))/eps);

  
  dcTdV(j) = 0.5*aLift*sigma*( eq6(mu2, lamNf2, thtaN) -  eq6(mu1, lamNf1, thtaN))/eps; % /sqrt(2);
  dcTdmu(j)= 0.5*aLift*sigma*( eq6(mu2, lamNf1, thtaN) -  eq6(mu1, lamNf1, thtaN))/eps; % /sqrt(2);
  

  
  xNM2=  ((eq9a(mu2, lamNf2, thtaN, oM(j))-eq11d(mu2, lamNf2, thtaN, oM(j)))-...
    (eq9a(mu1, lamNf1, thtaN, oM(j))-eq11d(mu1, lamNf1, thtaN, oM(j))))/eps;

  % the denominator on the right hand side...
  lamNf2=((V(j))*sin(alfaNf)-vi)/(oM(j)+eps)/rRotBs;
  mu2 = (V(j))*cos(alfaNf)/(oM(j)+eps)/rRotBs;
  % ....of the formula is the derivative with respect to oM
  xDN2=  ((eq9a(mu2, lamNf2, thtaN, oM(j))-eq11d(mu2, lamNf2, thtaN, oM(j)))-...
    (eq9a(mu1, lamNf1, thtaN, oM(j))-eq11d(mu1, lamNf1, thtaN, oM(j))))/eps;

  dcTdoM(j) = 0.5*aLift*sigma*( eq6(mu2, lamNf2, thtaN) -  eq6(mu1, lamNf1, thtaN))/eps;
  dcHdoM(j) = 0.5*sigma*aLift*(cH2sa(mu2,lamNf2,thtaN)- cH2sa(mu1,lamNf1,thtaN))/eps;
  

  doMdV2(j)= -xNM2/xDN2;

  % preparing the numerator... 
  lamNf3=((V(j))*sin(alfaNf+eps)-vi)/oM(j)/rRotBs;
  mu3 = (V(j))*cos(alfaNf+eps)/oM(j)/rRotBs;
  % ...for the alfa derivative
  xNM_alf=  ((eq9a(mu3, lamNf3, thtaN, oM(j))-eq11d(mu3, lamNf3, thtaN, oM(j)))-...
    (eq9a(mu1, lamNf1, thtaN, oM(j))-eq11d(mu1, lamNf1, thtaN, oM(j))))/eps;


  doMdalf(j)= -xNM_alf/xDN2;
  
  
  dTdV_Qn(j) = K_1*( dcTdV(j) + doMdV2(j)*(dcTdoM(j) + 2*cT/oM(j)));
  dTdV(j)    = K_1*( dcTdV(j));
  


  dAlf = eps/mu1 + 0.5*aLift*sigma*( eq6(mu3, lamNf3, thtaN)/sqrt(lamNf3^2 + mu1^2) -  eq6(mu1, lamNf1, thtaN)/sqrt(lamNf1^2 + mu1^2))/mu1;

  dcTdalf(j)=0.5*aLift*sigma*( eq6(mu3, lamNf3, thtaN)-eq6(mu1, lamNf1, thtaN))/dAlf;
  
  dTdAlf_Qn(j) = K_1*( dcTdalf(j) + doMdalf(j)*(dcTdoM(j) + 2*cT/oM(j)));
  dTdAlf_oM(j) = 8*K_1*( dcTdalf(j));
  
  % ======== dcTdV ============= 
  dmudV = mu/V(j);
  dmudalf = -mu*tan(alfaNf);
  dmudoM  = -mu/oM(j);
  
  
   
   
  dcTdlam = sigma*aLift/2*(0.5*bTl + 0.25*mu^2);
  dcTdmu  = sigma*aLift/2*(bTl*thtaN + 0.5*lamNf);
  
  
  m2l2= mu^2 + lamNf^2;
  s2m2l2 = 2*sqrt(m2l2);
  s3m2l2 = 2*m2l2^1.5;
  
  dlamdmu=((tan(alfaNf) - dcTdmu/s2m2l2 + mu*cT/s3m2l2)/...
           ( 1 +  dcTdlam/s2m2l2  - lamNf*cT/s3m2l2));
       
  dlamdmuAr(j)=dlamdmu;
  
  dlamdV = dlamdmu*dmudV;
    
  dcTdV(j) = dcTdmu*dmudV + dcTdlam*dlamdV;
  
  dTdV3_Qn(j) = K_1*( dcTdV(j) + doMdV2(j)*(dcTdoM(j) + 2*cT/oM(j)));
  dTdV3(j)    = K_1*dcTdV(j);

  % ======== dcTdAlf =========
  dlamdalf = dlamdmu*dmudalf;
  dlamdalfArry(j)=dlamdalf;
  
  dcTdAlf3(j) = dcTdmu*dmudalf + dcTdlam*dlamdalf;
  
  
  dcTdAlf2mu = 0.5*aLift*sigma*( eq6(mu1+eps, lamNf1, thtaN) -  eq6(mu1, lamNf1, thtaN))/eps;
  dcTdAlf2lam= 0.5*aLift*sigma*( eq6(mu1, lamNf1+eps, thtaN) -  eq6(mu1, lamNf1, thtaN))/eps;
  dcTdAlf2 = dcTdAlf2mu*dmudalf + dcTdAlf2lam*dlamdalf;

  dTdAlf2_Qn(j) = K_1*( dcTdAlf2 + doMdalf(j)*(dcTdoM(j) + 2*cT/oM(j)));
  dTdAlf2_oM(j) = K_1*dcTdAlf2;
  
  
  
  dTdAlf3_Qn(j) = K_1*( dcTdAlf3(j) + doMdalf(j)*(dcTdoM(j) + 2*cT/oM(j)));
  dTdAlf3_oM(j) = K_1*dcTdAlf3(j);
  
  cHmu(j) = 0.5*sigma*aLift*(cH2sa(mu+eps,lamNf,thtaN) - cH2sa(mu,lamNf,thtaN))/eps;
  cHlam(j)= 0.5*sigma*aLift*(cH2sa(mu,lamNf+eps,thtaN) - cH2sa(mu,lamNf,thtaN))/eps;
  
  
  dcHdV(j)= (cHmu(j)*dmudV + cHlam(j)*dlamdV);
  dHdV_oM(j)= K_1*(dcHdV(j) + doMdV(j)*(dcHdoM(j) + 2*cH(j)/oM(j)));

  dHdV_Qn(j) = K_1*dcHdV(j);
  
  % ========= dcHdalfa ================
  % Start of derivatives from naca-tn-D2655
  c1=sigma*aLift/8;
  c2=sqrt(mu^2+lamNf^2);
  cDNTR = (c2/oMR - vi*lamNf/c2/oMR^2 + c1/oMR);


  vHeli=mu*oMR/cos(alfaNf);
  
  c3 = (1-0.5*mu^2);
  c4 = (1+0.5*mu^2);

  

  dlamdalfaNf=8*mu^2*(bTl^2+1.5*mu^2)/((bTl^2-0.5*mu^2)*(8*mu+sigma*aLift));


  % dvi dthtaN
  dvidthtaN=c1*(2/3+mu^2)/cDNTR;
  dvidalfaNf=mu*(lamI^2/c2 + c1*(1-2*thtaN*mu*alfaNf))/cDNTR;

  da1dalfaNf = 2*mu/c3*(mu - dvidalfaNf/oMR) - a1 * c4/c3*tan(alfaNf);
  % dhcdmu=0.25*dDrgN*bTl^2;

  
  chi=lamNf/mu;
  
  xh= mu^2-mu*lamI*sin(alfaNf)/c2^2 + 2*mu^2/c3 - a1*tan(alfaNf)*c4/c3;
  xh= xh - mu/oMR*dvidalfaNf;
  dv1dalfaNf = dvidalfaNf*tan(0.5*chi) + vi/(1+cos(chi))*xh;

%   db1dthtaN = dv1dthtaN/oMR/c4;
%   fprintf('%s%10.5f\n','db1dthtaN ',da1dthtaN);

  db1dalfaNf =dv1dalfaNf/oMR/c3 - b1*mu^2* c4/c3*tan(alfaNf);
  %fprintf('%s%10.5f\n','db1dalfNf ',da1dalfaNf);

%   db1dmu = (dv1dV-b1Br*mu*cos(alfaNf))/c3;
%   fprintf('%s%10.5f\n','db1dmu    ',db1dmu);

  
  dcHdalfaNf = mu*(thtaN*lamNf - dDrgN/aLift - 0.5*a1^2 -b1*lamI/8)*tan(alfaNf) - thtaN*mu^2 +1.5*a1*mu;
  dcHdalfaNf = dcHdalfaNf + dvidalfaNf/oMR*(thtaN*mu - 1.5*a1);
  dHdalfaNf(j) = K_1*c1*(dcHdalfaNf + da1dalfaNf*(2/3*thtaN + 1.5*lamNf + a1*mu) + db1dalfaNf*mu*lamI/8 + dv1dalfaNf*b1*mu/oMR/8);
  % fprintf('%s%10.5f%s%10.5f\n','dcHdalfaNf ',dcHdalfaNf/sigma,'   dcHdalfaNf ',dcHdalfaNf);

  
  
  dcHdalf(j)= -(cHmu(j)*dmudalf + cHlam(j)*dlamdalf);
  dHdalf_oM(j) = K_1*dcHdalf(j);
  dHdalf_Qn(j) = K_1*(dcHdalf(j) + doMdalf(j)*(dcHdoM(j) + 2*cH(j)/oM(j)));
  
  % eq40
  dlta= -5.5e-4*V(j) + 0.042;
  
  
  dMtdV_f12(j)    = 1.9*mu - 0.5;
  dMnpdalf_f12(j) =-390*mu + 16;
  dMtdalf_f12(j)  =-1063.8*mu + 77.66;
  
  
  alf_t=alfaNf -dlta -iRigMR;
  alf_p=alf_t;
  
  M_r_alf = (dMnpdalf_f12(j) + dMtdalf_f12(j))*alf_t;
  
  
  a1s = eq2(mu, lamNf, thtaN) - dlta;

  % rotor force leaver arms
  alfaS=(alfaD + a1s);
  lOne=  xCG          - hMR*alfaS;
  hOne= -xCG*alfaS + hMR;
  cMs = 1.0e-4;
  
  
  lamD=mu*alfaD - lamI;
  ca1Br = (1+0.5*mu);
  
  da1dmu=2*(4/3*bTl*thtaN+lamD)*(bTl^2 - 1.5*mu^2)/(bTl^2 +1.5*mu^2)^2*ca1Br;
  
  mU1 = -hOne*(dHdV_oM(j) + mHeliG*alfaS) + lOne*dTdV3_Qn(j) + cMs*da1dmu;

  Mu(j) = mU1 + dMtdV_f12(j);
%  fprintf('%s%10.6f%s%10.6f%s%10.6f\n','mU     ',mU1,' + M_r_a ',M_r_alf,'  = ',Mu(j));

%   mW1 = -hOne*xW +lOne*zW + cMs*da1dw + mwf;
%   mW  = mW1+ mWT;
%   fprintf('%s%10.6f%s%10.6f%s%10.6f\n','mW     ',mW1,' + mWT ',mWT,'  = ',mW);
% 
%   mQ1 = -hOne*xQ + lOne*zQ - 16*cMs/(gamRot*(1-0.5*mu^2)) + mqf;
%   mQ1 = -hOne*xQ + lOne*zQ +cMs*da1dq + mqf;
% 
%   mQ = mQ1 +mQT;
%   fprintf('%s%10.6f%s%10.6f%s%10.6f\n','mQ     ',mQ1,' + mQT ',mQT,'  = ',mQ);

  
end


figure; plot(muAr, lamAr); title('lam vs mu'); grid

oM= oM/2/pi*60;
figure; plot(muAr,oM); title('oM vs mu'); grid
h=get(0,'CurrentFigure');
saveas(h,'test1','tif');

figure; plot(muAr,3.6*V/1.605); title('V[mph] vs mu')
figure; plot(V/f2m,rad2deg(alNf)); title('alNf vs V[ft/s]'); grid
xh= 1.0; % 0.95/1.27;
figure; plot(muAr,xh*doMdV*f2m); hold on
plot(muAr,xh*doMdV2*f2m,'g'); title('fig. 5  doM/dV vs mu'); legend('anly','nmrcl');grid

figure; plot(muAr,doMdalf); title('fig. 5  doM/dalf vs mu');grid
cXf  = f2m/lbForce;  % *18/24.5;
figure; plot(muAr,dTdV_Qn*cXf);   hold on; plot(muAr,dTdV*cXf,'g');  title('dTdV vs mu');  legend('Q=0','oM=C'); grid
figure; plot(muAr,dTdV3_Qn*cXf);  hold on; plot(muAr,dTdV3*cXf,'g'); title('dTdV3/Qn vs mu');  legend('Q=0','oM=C'); grid


% dTdAlf;
cX  = 1/lbForce;

dTdAlf_Qn=dTdAlf_Qn*cX; dTdAlf_oM=dTdAlf_oM*cX;
dTdAlf2_Qn=dTdAlf2_Qn*cX; dTdAlf2_oM=dTdAlf2_oM*cX;
dTdAlf3_Qn=dTdAlf3_Qn*cX; dTdAlf3_oM=dTdAlf3_oM*cX;


figure; plot(muAr,dTdAlf3_Qn); hold on; plot(muAr,dTdAlf3_oM,'g'); title('fig. 9 dTdAlf3 vs mu'); legend('Q=0','oM=C');grid 
figure; plot(muAr,dTdAlf_Qn); hold on; plot(muAr,dTdAlf_oM,'g'); title('fig. 10 dTdAlf vs mu'); legend('Q=0','oM=C');grid
figure; plot(muAr,dTdAlf2_Qn); hold on; plot(muAr,dTdAlf2_oM,'g'); title('fig. 10a dTdAlf2 vs mu'); legend('Q=0','oM=C'); grid

% H-force
figure; plot(muAr,K_1*cH*cX); title('fig. 8 H-force lb');grid
figure; plot(muAr,dHdV_Qn*cXf);  hold on; plot(muAr,dHdV_oM*cXf,'g');   
title('fig. 10 dH/dV - oMlb/(ft/sec)'); grid; legend('dH/dV Qn','dH/dV oM')

figure; plot(muAr,dHdalf_Qn*cX); hold on; plot(muAr,dHdalf_oM*cX,'g');  plot(muAr,dHdalf_26*cX,'y');
legend('dH/dalf Qn','dH/dalf oM'); title('fig. 11 dH/dalf');grid

figure; plot(muAr,dHdalfaNf);
legend('dH/dalf 2665'); title('fig. 11 dH/dalf  2665');grid

figure; plot(muAr,dMtdV_f12*100); title('fig. 11 dM r t p '); hold on
plot(muAr,dMnpdalf_f12,'g');   plot(muAr,dMtdalf_f12,'r'); grid; legend('dMtdV','dMnpdalfx100','dMtdalf');

figure; plot(muAr,Mu); title('fig. 12 dM/dV');grid

figure; plot(muAr,cTsAr); title('fig. 13 cT vs mu');grid

figure; plot(muAr,dcTdAlf3/sigma); title('fig. 14 dcTdAlf3 vs mu');grid
figure; plot(muAr,alfaDegAr); title('fig. 15 alfaDeg vs mu');grid

figure; plot(muAr,lam_Hov); title('fig. 16 lam_Hov vs mu');grid
figure; plot(mu_lamHov,lam_Hov); title('fig. 17 lam_Hov vs mu_lam');grid


% dcHsalfaNf    0.01005   dcHdalfaNf    0.00036
% http://www.aero.psu.edu/avia/pubs/YomHor12.pdf
% http://www.icas.org/ICAS_ARCHIVE/ICAS2006/PAPERS/017.PDF
% http://www.youtube.com/watch?v=iRT62HT3aiM
% http://www.ebay.com/itm/1940-Harold-Pitcairn-Autogiro-Co-Works-on-Rotor-Hub-of-PA-36-Wire-Photo-/300777663002?afsrc=1
% http://justacarguy.blogspot.de/2011/06/want-to-draw-crowd-put-helicopter-in.html
% http://www.rotaryforum.com/forum/showpost.php?p=344834&postcount=1   Chuck on Profiles
% http://www.mashpedia.com/Pitcairn_PCA-2?tab=1&startvid=40
% http://www.whatifmodelers.com/index.php?action=printpage;topic=32392.0
% http://www.reaa.ru/cgi-bin/yabb/YaBB.pl?num=1264932310
% http://forum.woodenboat.com/showthread.php?80245-Tall-Ships/page1
% http://www.gutenberg.org/files/44471/44471-h/44471-h.htm  BristolPrivateer
% http://themodelshipwright.blogspot.de/2013/02/free-ship-plan-algerian-xebec-two-fer.html
% http://www.britishpathe.com/video/windmill-girls-at-seaside-issue-title-pathe-pictor
% http://www.britishpathe.com/gallery/queen-surprising-facts/4
% http://www.britishpathe.com/video/speed-and-comfort
% http://www.mashpedia.com/Pitcairn_Aircraft_Company?pagetype=search&tab=1&startvid=47&pagecode=CGQQAA&xn=101&autoplay=1 Air Mail R-4B
% http://www.mashpedia.com/Pitcairn_Aircraft_Company?pagetype=search&tab=1&startvid=8&pagecode=CGQQAA&xn=101&autoplay=1
% bensen towed gyro 
% http://www.britishpathe.com/video/man-makes-his-own-autogiro/query/prescott
% cross country Martin B-10
% http://www.mashpedia.com/Pitcairn_Aircraft_Company?pagetype=search&tab=1&startvid=46&pagecode=CGQQAA&xn=101&autoplay=1
% http://www.youtube.com/watch?v=Mjz8vwDUQpk  Battlefield Atlantic
% http://www.youtube.com/watch?v=nPxZiZVRBhE  Falklands
% http://www.youtube.com/watch?v=FKvRHaJ2w6w  Ufberth
% Lockheed 5B (5C) Vega NR926Y c/n 134 "Lituanica II" (ex Shell Oil Co.), a specially modified Vega used for a non-stop flight to Lithuania from New York on September 21-22, 1935, flown by pilot Lt. Felix Waitkus. Needless to say it was filled with fuel tanks, and the back windows were blanked over. The airplane suffered a forced landing at Ballinrobe, County May, Ireland on the 22nd. The right landing gear collapsed and the wing and fuselage were damaged. Pilot Waitkus was uninjured. It was rebuilt in Lithuania, to Lithuanian Air Corps. 
% http://www.flickr.com/photos/kemon01/11969520834/in/photostream/
% https://www.youtube.com/watch?v=hUsOQyfnfDY Air Battle Guadalcanal


