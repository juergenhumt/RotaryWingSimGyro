function pwrOut = pwrCalc(height, mu, oM, trPwr, pwrCnstHeight, nOut)
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
%         --- Version 1.2 ---
%
% input:
%   height     height at which aircraft is flying
%   mu         tip speed ratio
%   oM         rotational speed of rotor
%   trPwr      ratio of available to power consumed by tail rotor
%   pwrCnstHeight  height up to which engine power is constant
%   nOut       ratio of hover hight to rotor radius
%   
% output:
%   pwrOut     power absorbed by main rotor  
% 
% H  = [  0 m 1000 m 2000 m 3000 m 4000 m 4500 m 5000 m ]
% cL = [ 0.654 0.723 0.799 0.885 0.981 1.034 1.090 ]
% cD = [   0.0120 0.0126 0.0134 0.0147 0.0181 0.0250 0.0393 ]

global sigma aDiskMR rhoAir rRot thtaClmb fpA...
    mHeliG engineRPM engineHP cLpwrCff cDpwrCff f2m

% calculate air density at given   height 
[t, p, rhoAirH, a] = atmosphere(height);



rhoAirN= 1.2252;

oMR = oM*rRot;
vHn = mu*oMR;
mu2 = mu*mu;

cT = mHeliG/(rhoAirH*aDiskMR*oMR^2);

% climbing speed is caculated from climbing angle
wClmb = thtaClmb*vHn;

% power available in dimensions of Watt [W]
pAV = engineHP*0.743*1000;

% engine power may be constant up to some height
if height > pwrCnstHeight
  pAV = pAV*(1.11*rhoAirH/rhoAirN - 0.11);
else 
  if height > 1.0
    pAV = pAV*(1.013*rhoAirH/rhoAirN - 0.013);
  end
end
cPav = pAV/(rhoAirH*aDiskMR*oMR^3);

lamHvr = sqrt(0.5*cT);

% nOut is used to give ratio of hover hight to rotor radius
% for hovering rotor
if (mu < 1.0e-3) & (nOut > 0.55)
   lamHvr = lamHvr*(1 - 0.5/(1 + 4.0*(1/nOut)^2)); 
end    



vHvr = lamHvr*oMR;


p = (vHn/vHvr)^2;

cL = 6*cT/sigma;

% xh  = crvVal(thtaN, cLpwrCff);
cDN = crvVal(cL, cDpwrCff);


muX=0.0005;
% induced power
if mu < muX
  kP   = 1.15 + 0.05*mu/muX;
  lamI = lamHvr*(1 - mu/muX  + mu/muX*lam4slv(p));
else
  kP   = 1.20;
  lamI = lamHvr*lam4slv(p);
end

% Westland Helicopter ltd constant, similar
% to the one used by Bramwell
kWstlnd = 4.65; 

rDrag = 0.5*rhoAirH*fpA*vHn*vHn;

zeta  = 1.3;
ztaDsc= 1.0;

cPreq = kP*lamI*cT + 0.125*sigma*cDN*(1 + kWstlnd*mu2) + mu*rDrag/mHeliG*cT + zeta*wClmb/oMR*cT;

pReq = (1.0+trPwr)*cPreq*rhoAirH*aDiskMR*oMR^3;

pReqHp = pReq/1000/0.743;

zNN = (zeta*cT);

lamDsc  = - kP*lamI/ztaDsc - 0.125*sigma*cDN/(ztaDsc*cT)- 0.5*mu*mu2*fpA/((ztaDsc*cT)*aDiskMR);

lamClmb = (1.0 - trPwr)*cPav/zNN - kP*lamI/zeta - 0.125*sigma*cDN/zNN - 0.5*mu*mu2*fpA/(zNN*aDiskMR);

vDsc  = lamDsc*oMR;
vClmb = lamClmb*oMR;

if nOut < 0.5
  pwrOut = [vClmb, vDsc, vHn, cT, pReq];
else
  pwrOut = abs(vClmb);
end
 
