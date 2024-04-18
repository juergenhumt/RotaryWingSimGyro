function [aN, a1, a2, b1, b2, dQ] = clcRotAnglTwst(mu, lamNf, thtaN, oMR, rotThrst)
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
%         --- Version 2.0 ---
%
%
% This module calculates twist coefficients for rotor blades with offset between elastic
% and aerodynamic axis using formulae in naca 591 for a give flight state 
%
% Input:
%   mu      = advance ratio  -> V*cos(alfaD)/oMR
%   lamNf   = inflow at no feathering axis (control/swash plate)
%   thtaN   = collective pitch at blade root
%   oMR     = tip speed ratio
%   rotTrst = thrust currently developed by rotor
%
% Output:
%   a0      = rotor coning angle
%   a1      = first harmonic longitudinal coefficient
%   a2      = second harmonic longitudinal coefficient
%   b1      = first harmonic lateral coefficient
%   b2      = second harmonic lateral coefficient 
%   dQ      = difference between de- and accellerating torque
%

global gamRot bTl thta1 aLift dDrgN rRot cMIN


% calculate twist infulence coefficients
[eNTwst, e1Twst, e2Twst, n1Twst, n2Twst] = clcCoeffTwst(mu, lamNf, thtaN, oMR, rotThrst);
xNNab=(gamRot^2*bTl^8+144);

a2n= mu^2*gamRot/xNNab*(lamNf*bTl*(16 + 7/108*gamRot^2*bTl^8) + thtaN*bTl^2*(46/3 + 7*gamRot^2*bTl^8/144) +...
    thta1*bTl^3*(12 + 7*gamRot^2*bTl^8/180));
    
da2 = mu^2*gamRot/xNNab*(eNTwst*bTl^3*(12 + 7*gamRot^2*bTl^8/180) - e1Twst/mu/30*gamRot*bTl^8 + 2/5*n1Twst/mu*bTl^4 +...
    24/5*e2Twst/mu^2*bTl*5 + 2/5*n2Twst/mu^2*gamRot*bTl^9);

a2 = a2n + da2;
% b2 = -mu^2*gamRot^2/(gamRot^2*bTl^3 + 144)*(5/9*lamNf*bTl^5 + 25/36*thtaN*bTl^6 + 8/15*(thta1 + eNTwst)*bTl^7  +...
%      2/5/gamRot*e1Twst/mu*bTl^4 + n1Twst/mu/30*bTl^8 + 2/5*e2Twst/mu^2*bTl^9 - 24/5/gamRot*n2Twst/mu^2*bTl^5);


b2n = 5/9*lamNf*bTl^5 + 25/36*thtaN*bTl^6 + 8/15*thta1*bTl^7;
db2 = 8/15*eNTwst*bTl^7  + 2/5/gamRot*e1Twst/mu*bTl^4 + n1Twst/mu/30*bTl^8 + 2/5*e2Twst/mu^2*bTl^9 - 24/5/gamRot*n2Twst/mu^2*bTl^5;

b2 = -mu^2*gamRot^2/xNNab*(b2n + db2);

% aN ------ aN
aNn = 0.5*gamRot*(lamNf*(bTl^3/3 + 0.080*mu^3) + thtaN*(0.25*bTl^4 + 0.25*mu^2*bTl^2 - mu^4/32) +...
        thta1*(0.2*bTl^5 + mu^2*bTl^3/6) + 0.125*mu^2*b2*bTl^2) - cMIN/oMR^2; % cMI= Mblade/(I1blade*(oM*rRot)^2)
    
daNn = 0.5*gamRot*(eNTwst*(0.2*bTl^5 + mu^2*bTl^3/6) + 0.25*mu*n1Twst*bTl^4 - mu^2*e2Twst*bTl^3/12);

aN = aNn + daNn;

xNNa=bTl^4 - 0.5*mu^2*bTl^2;


% a1 ------ a1
a1n = 2*mu/xNNa*(lamNf*(bTl^2-0.25*mu^2) + thtaN*(4/3*bTl^3 + 0.106*mu^3) + thta1*bTl^4 -...
    b2*bTl^3/3);

da1= 2*mu/xNNa*(eNTwst*bTl^4 + n1Twst/mu*(2/5*bTl^5 + 0.5*mu^2*bTl^3) -0.5*e2Twst*bTl^4);

a1 = a1n + da1;

% b1 ------ b1
xNNb=bTl^4+0.5*mu^2*bTl^2;

b1n= 4*mu/xNNb*(aN*(bTl^3/3 + 0.35*mu^3) + a2*bTl^3/6);

db1 = 4*mu/xNNb*(- e1Twst*(0.2*bTl^5+mu^2*bTl^3/12) - 0.25*n2Twst*bTl^4);

b1 = b1n + db1;


aF = 0.5*bTl^2 - 0.25*mu^2;

bF = thtaN*bTl^3/3 + 2/9/pi*mu^3*thtaN + 0.25*thta1*bTl^4 + mu^4*thta1/32 +...
     mu*a1*(0.5*bTl^2 - 3/8*mu^2) + eNTwst*(0.25*bTl^4 + mu^4/32) +...
     0.5*n1Twst*mu*bTl^3/3;


cF = aN^2*(0.25*mu^2*bTl^2 - mu^4/16) - mu*aN*b1*bTl^3/3 +...
     a1^2*(0.125*bTl^4 + 3/16*mu^2*bTl^2) + b1^2*(0.125*bTl^4 + mu^2*bTl^2/16) -...
     a2*(0.25*mu^2*aN*bTl^2 + mu*b1*bTl^3/6) + 0.5*a2^2*bTl^4 +...
     b2*(0.125*mu^2*thtaN*bTl^2 + mu^2*thta1*bTl^3/12 + mu*a1*bTl^3/6) + 0.5*b2^2*bTl^4 -...
     0.25*dDrgN/aLift*(1.0 + mu^2 - 0.125*mu^4) + eNTwst*(mu^2*b2*bTl^3/12)+...
     0.5*e1Twst*(-0.25*mu*aN*bTl^4 + b1*(0.2*bTl^5 + mu^2*bTl^3/12) - 0.125*mu*a2*bTl^4) -...
     0.5*n1Twst*(-a1*(0.2*bTl^5 - mu^2*bTl^3/12) - 0.125*mu*b2*bTl^4) +...
     0.5*e2Twst*(0.25*mu*a1*bTl^4 + 2/5*b2*bTl^5) +...
     0.5*n2Twst*(-mu^2*aN*bTl^3/6 + 0.25*mu*b1*bTl^4 - 2/5*a2*bTl^5);
 
 
 % cF = cF + 1.2e-2*mu*mu - 3.7e-4;
 % cF = cF + 3.1e-2*(0.02-lamNf);
 cF = cF + 1.10*lamNf*(0.037-lamNf);

 
 
 dQ = aF*lamNf^2 + bF*lamNf + cF;

 
