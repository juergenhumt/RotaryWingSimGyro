function [alfaPr, cTs, DvLN] = PvLMuIt(mu, thtaN, vErth, PvL, oM)
global bTl gamRot rhoAir aLift thta1 cMI  sigma fpA rhoAir vHeli mHeliG...
    alfaNf oM rRot muGlo PvLGlo cDrgGlo thtaClmb...
%     t11 t12 t13 t14 t15 t16 t17 t18 t19 t11N ...
%        t21 t22 t23 t24 t25 t26 t31 t32 t33 dt31dmu dt32dmu dt33dmu...
%        t41 t42 t43 t44 t45 t46 t47 t48 t49 t41N...
%        t51 t52 t53 t54 t55 t56 t57 t58 t59 t51N t511 t512 t513 t514...
%        t61 t62 t63 t64 t65 t66 t67 t68 t69 t61N t611 t612 t613 t614...    
%        t71 t72 t73;
%
%
% this module calculates the angle a prime (alfaPr) as given by equation (5)
% page 14 % derived on pages 13/14 in naca 2309. The definition of a prime 
% in the report is:
% projection of angle between rotor force vector and axis of no feathering
% in plane containing flight path and angle of no feathering
%
%

warning off
% reynolds number correction
if (PvL < 0)
   xS = [0.25];
else
   xS = [-0.01];
end 

% initial value
xS= -0.05;
vNorm = norm(vErth);

options = optimset('MaxFunEvals',500,'TolFun',1.0e-9); % ,'Display','off');
% cQit(lamNf, mu, thtaN, vGyro, DvLp, PvL)
DvLp = 0.5*rhoAir*norm(vErth)^2*fpA/mHeliG;
DvLc = tan(thtaClmb);

% Matlab
[lAmNf,fVal] = fsolve(@cQit,xS,options,mu,thtaN, DvLp, DvLc, PvL, oM);
% octave


% rotor coefficients
aN = eq1(mu,lamNf,thtaN);
a1 = eq2(mu,lamNf,thtaN);
b1 = eq3(mu,lamNf,thtaN);

% a2 = mu^2*(t21*lamNf+t22*thtaN+t23*thta1);
% b2 = mu^2*(t24*lamNf+t25*thtaN+t26*thta1);
% 
lamP = lamNf+mu*a1;

% equation 6 page 212
% 2*cT/(sigma*aLift) = t31...
cT2sA   = eq6(mu, lamNf, thtaN);

% equation on lower left hand side of page 218
cQ2sA  = eq9a(lamNf,thtaN); % fprintf('%s%10.4f\n','cQ2sA      ',cQ2sA);

cQ2sd = eq11d(mu, lamNf,thtaN);

cTs = 0.5*cT2sA*aLift;
cT  = cTs*sigma;

% equation 14 page 215 mu*2*cT/(Sigma*aLift)=...
mucT2sA = eq14(lamNf, thtaN);

% Rotor drag lift ratio. 
% Equation on upper right of page 219
DvLN = eq13(lamNf,thtaN)/mucT2sA;

cNN = (2*mu*(mu^2+lamNf^2)^0.5);

% fprintf('%s%10.4f%s%10.4f\n','lAmIT      ',xOut,'   lamNf  ',lamNf);

% DvL2 = sigma*dDrgN*(1+3*mu^2+3/8*mu^4)/8/mu/cT; % +cT/cNN;

DvLi = cT/cNN;

DvLp = fpA*rhoAir/2*(vHeli)^2/mHeliG;
DvLc = atan(thtaClmb);

alfaNf = atan(lamP/mu+DvLi);

% equation 26 page 218
cLs = cT2sA*aLift*cos(alfaNf)^3/mu^2;

DvLu = DvLp + DvLu;

 % eq. (A10) p. 32 NACA2309
alfaPr = 57.3*atan((DvLN-PvL+(t32*thtaN-cT2sA)/mu/t31)/(1-DvLu*tan(alfaNf)));
