function outRes = heliDeriv2655(mu, alfaNf, lamNf, lamI, thtaN, a1, cT, cH, oMR, vHeli)
global bTl rRot sigma thta1 gamRot aLift dDrgN dDrg1 dDrg2



nU=0.7;

vi=lamI*oMR;
aPr=atan(cH/cT);

aN = eq1(mu,lamNf,thtaN);


ch1=atan(real(-mu/lamNf));
if ch1 > 0.5*pi
  ch1 = pi-ch1;
end


chi = abs(ch1+a1);
chi2= 0.5*chi;

v1= vi*tan(chi2);
lamPr = v1/oMR;

% Start of derivatives from naca-tn-D2655
c1=sigma*aLift/8;
c2=(mu^2+lamNf^2);
c2s=sqrt(c2);
cDNTR = (c2s/oMR - vi*lamNf/c2s/oMR^2 + c1/oMR);


% dvi dthtaN
dvidthtaN=c1*(2/3+mu^2)/cDNTR;
% fprintf('%s%10.5f\n','dvidthtaN ',dvidthtaN);

dvidalfaNf=mu*(vi^2/c2s/oMR^2 + c1*(1-2*thtaN*mu*tan(alfaNf)))/cDNTR;
% fprintf('%s%10.5f\n','dvidalfaNf ',dvidalfaNf);

% dcTsdmu = 2*aLift*mu/(8*mu + aLift*sigma)*(2*thtaN*mu*cos(alfaNf) + sin(alfaNf) + 0.5*cT/mu^2);
dvidV = (c1*(2*thtaN*mu*cos(alfaNf)+sin(alfaNf)) - vi*(vHeli - vi*sin(alfaNf))/c2s/oMR^2)/(cDNTR*oMR);
% fprintf('%s%10.5f\n','dvidV ',dvidV);

cx1 = 1 + cos(chi);
c3=1-0.5*mu^2;
c4=1+0.5*mu^2;

% b1 = (mu*aN + 1.1*nU^0.5*lamI)/c3;
b1 = eq3(mu, lamNf, thtaN);

dv1dthtaN= dvidthtaN*tan(chi2) + vi/cx1*(8/3*mu/c3 -(1/c2 + 2/c3)*mu/oMR*dvidthtaN);

dv1dalfaNf= dvidalfaNf*tan(chi2) + vi/cx1*((vHeli^2- vHeli*vi*sin(alfaNf))/oMR^2/c2 +2*mu^2/c2 - c4/c3*a1*tan(alfaNf)...
    -(1/c2 + 2/c3)*mu/oMR*dvidalfaNf);
   

dv1dV=dvidV*tan(chi2) + vi/cx1*((vi*cos(alfaNf))/oMR^2/c2 + c4/c3*a1*cos(alfaNf)/mu/oMR...
    + 2*mu*sin(alfaNf)/c3/oMR - (1/c2 + 2/c3)*mu/oMR*dvidV);



da1dthtaN= 2*mu/oMR*(4/3*oMR - dvidthtaN)/c3;

da1dalfaNf= 2*mu/c3*(mu - dvidalfaNf/oMR) - a1*c4/c3*tan(alfaNf);

da1dmu = a1*c4/mu/c3*cos(alfaNf) + 2*mu/c3*(sin(alfaNf) - dvidV);

db1dthtaN=dv1dthtaN/oMR/c4;

db1dalfaNf=dv1dalfaNf/oMR/c4  + b1*mu^2*tan(alfaNf)/c4;

db1dmu = (dv1dV - b1*mu*cos(alfaNf))/c4;

dcTdthtaN= 0.25*sigma*aLift*(2/3 + mu^2 - dvidthtaN/oMR);

dcTdalfaNf= 0.25*sigma*aLift*(mu + 2*thtaN*mu^2*tan(alfaNf) - dvidalfaNf/oMR);

% dcTdmu = 2*aLift*sigma*mu/(8*mu + aLift*sigma)*(2*thtaN*mu*cos(alfaNf) + sin(alfaNf) + 0.5*cT/mu^2);
% dvidV = (c1*(2*thtaN*mu*cos(alfaNf)+sin(alfaNf)) - vi*(vHeli - vi*sin(alfaNf))/c2s/oMR^2)/(cDNTR*oMR);
dcTdmu= 0.25*sigma*aLift*2*thtaN*mu*cos(alfaNf) + sin(alfaNf) + dvidV;


dcHdthtaN= 0.25*sigma*aLift*(2/3*a1 - mu*lamNf - dvidthtaN*(thtaN*mu - 1.5*a1)/oMR + da1dthtaN*(1.5*lamNf + a1*mu+ 2/3*thtaN)+...
    db1dthtaN*mu*lamPr/8 + dv1dthtaN*b1*mu/8/oMR);

dcHdalfaNf= 0.25*sigma*aLift*(mu*(thtaN*lamNf -dDrgN/aLift -0.5*a1^2 -b1*lamPr/8)*tan(alfaNf)...
    - thtaN*mu^2 + 1.5*a1*mu + (thtaN*mu -1.5*a1)*dvidalfaNf/oMR...
    +da1dalfaNf*(2/3*thtaN+1.5*lamNf + a1*mu) + db1dalfaNf*mu*lamPr/8 +dv1dalfaNf*b1*mu/oMR/8);

dcHdmu= 0.25*sigma*aLift*(dDrgN*cos(alfaNf)/aLift -thtaN*(lamNf*cos(alfaNf) + mu*sin(alfaNf)) + 1.5*a1*sin(alfaNf)...
        +0.5*a1^2*cos(alfaNf) + b1*lamPr*cos(alfaNf)/8 - dvidV*(thtaN*mu -1.5*a1)...
        +da1dmu*(2/3*thtaN +1.5*lamNf+a1*mu) + db1dmu*mu*lamPr/8 + dv1dV*b1*mu/8);

    
dcQdthtaN= 0.5*sigma*aLift*(dv1dthtaN/oMR*0.25*(b1 -lamPr) - db1dthtaN*(0.25*(b1-lamPr) + b1*mu^2/8)...
    - lamNf/3 + dvidthtaN/oMR*(thtaN/3 + lamNf + 0.25*a1*mu) - da1dthtaN*(0.5*mu*lamNf + 0.25*a1 + 3/8*a1*mu^2));

dcQdalfaNf= 0.5*sigma*aLift*(mu*(vi/oMR-thtaN/3) - mu^2*(0.5*dDrgN/aLift + 1.0 - 3/8*a1^2 - b1^2/8 -0.25*a1*tan(alfaNf))*tan(alfaNf)...
     - 0.5*a1*mu^2 - 0.5*a1*vi*mu/oMR*tan(alfaNf) + dvidalfaNf/oMR*(thtaN/3 +lamNf + 0.5*a1*mu) + da1dalfaNf*(0.5*mu*vi/oMR - a1*(0.25 + 3/8*mu^2)...
     - 0.5*mu^2*tan(alfaNf)) + db1dalfaNf*(0.25*lamPr - b1*(0.25 + mu^2/8)) + dv1dalfaNf*0.25*(b1-lamNf));
 
dcQdmu = 0.5*sigma*aLift*(0.5*mu*(dDrgN/aLift - 3/4*a1^2 - 0.25*b1^2)*cos(alfaNf) - (thtaN/3 +lamNf + a1*mu)*sin(alfaNf)...
     + 0.5*a1*vi/oMR*cos(alfaNf) + dvidV*(thtaN/3 + lamNf + 0.5*a1*mu) - da1dmu*(a1*(0.25 + 3/8*mu^2) + 0.5*mu*lamNf)...
     + db1dmu*(0.25*lamPr - 0.25*b1*(1+0.5*mu^2)) + dv1dV*0.25*(b1-lamPr));

daPrdthtaN= (dcHdthtaN - dcTdthtaN*tan(aPr))*cos(aPr)^2/cT;
daPrdalfaNf= (dcHdalfaNf - dcTdalfaNf*tan(aPr))*cos(aPr)^2/cT;
daPrdmu    = (dcHdmu- dcTdmu*tan(aPr))*cos(aPr)^2/cT;
 
outRes =[ aN a1 b1 chi real(da1dthtaN) real(da1dalfaNf) real(da1dmu)....
    real(db1dthtaN) real(db1dalfaNf) real(db1dmu)...
    real(dcTdthtaN) real(dcTdalfaNf) real(dcTdmu)...
    dcHdthtaN dcHdalfaNf dcHdmu...
    dcQdthtaN dcQdalfaNf dcQdmu...
    daPrdthtaN daPrdalfaNf daPrdmu dvidV];