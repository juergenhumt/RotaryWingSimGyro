function [eNTwst, e1Twst, e2Twst, n1Twst, n2Twst] = clcCoeffTwst(mu, lamNf, thtaN, oMR, rotThrst)
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
% This module calculates twist coefficients for rotor blades with offset
%  between elasticand aerodynamic axis formulae from naca-tn-600
%
% Input:
%   mu      = advance ratio  -> V*cos(alfaD)/oMR
%   lamNf   = inflow at no feathering axis (control/swash plate)
%   thtaN   = collective pitch at blade root
%   oMR     = tip speed ratio
%   rotTrst = thrust currently developed by rotor
%
% Output:
%   eNTwst, e1Twst, e2Twst, n1Twst, n2Twst: all twist coefficients used in naca 600
%



global bTl aLift gamRot rRot nB thta1 cLftGlo cDrgGlo sigma...
  thta1 aLift f2m Gmod...
  lbForce lbMass rhoAir gEarth cGA CmPrfl cHMR...
  
cRed = 1.0;

% eNTwst = rad2deg(rotThrst*cGA/nB/Gmod + 0.5*rhoAir*oMR^2*rRot*cHMR^2*CmPrfl*(bTl^3/3 + 0.5*mu^2*bTl)/Gmod);

Axx = 0.5*rhoAir*cHMR*aLift*oMR^2*rRot*cGA/Gmod;


c_eN1= 0.5;  % 0.5 ++!!!++
eNTwst = rotThrst*cGA/nB/Gmod;
eNTwst2= Axx*CmPrfl*cHMR/aLift/cGA*(bTl^3/3 + c_eN1*mu^2*bTl);

eNTwst = (eNTwst + eNTwst2);

e1Twst = -mu*gamRot*Axx*(lamNf*(bTl^5/108 - 0.0161*mu^2*bTl^3) + thtaN*(bTl^6/144 - 0.0049*mu^2*bTl^4) +...
     (thta1+eNTwst)*(bTl^7/180 - 0.0048*mu^2*bTl^5)) +...
     mu^3*gamRot*Axx^2*(0.0021*lamNf*bTl^7 + 0.0007*thtaN*bTl^8 - 0.0002*(thta1 + eNTwst)*bTl^9 +...
     0.0071*CmPrfl*cHMR/aLift/cGA*bTl^8);
 
n1Twst = mu*Axx*(lamNf*(bTl/3 + 0.341*mu^2/bTl) + thtaN*(bTl^2/9 + 0.233*mu^2) + 0.175*mu^2*(thta1 + eNTwst)*bTl + CmPrfl*cHMR/aLift/cGA*bTl^2) +...
     mu^3*Axx^2*(0.007*lamNf*bTl^3 + 0.006*thtaN*bTl^4 + 0.005*(thta1+eNTwst)*bTl^5 - 0.003*CmPrfl*cHMR/aLift/cGA*bTl^4);

e2Twst = mu^2*Axx*(0.796*lamNf + 0.578*thtaN*bTl + 0.554*(thta1+eNTwst)*bTl^2 - 0.5*CmPrfl*cHMR/aLift/cGA*bTl) -...
      mu^2*Axx^2*(0.032*lamNf*bTl^4 + 0.044*thtaN*bTl^5 + 0.039*(thta1+eNTwst)*bTl^6 - 0.055*CmPrfl*cHMR/aLift/cGA*bTl^5);
  
  
n2Twst = -mu^2*gamRot*Axx*(0.0291*lamNf*bTl^4 + 0.0290*thtaN*bTl^5 + 0.0225*(thta1+eNTwst)*bTl^6) -...
      mu^2*gamRot*Axx^2*(0.0078*lamNf*bTl^8 + 0.0057*thtaN*bTl^9 + 0.0052*(thta1+eNTwst)*bTl^10 - 0.0050*CmPrfl*cHMR/aLift/cGA*bTl^9);
  
  
% e3 =  mu^3*gamRot*Axx*(0.0085*lamNf*bTl^3 + 0.0143*thtaN*bTl^4 + 0.0108*(thta1+eNTwst)*bTl^6) +...
%       mu^3*gamRot*Axx^2*(0.0135*lamNf*bTl^7 + 0.0111*thtaN*bTl^8 + 0.0097*(thta1+eNTwst)*bTl^9 - 0.0047*CmPrfl*cHMR/aLift/cGA*bTl^8);
% 
% n3 =  mu^3*Axx*(0.271*lamNf/bTl + 0.381*thtaN + 0.280*(thta1+eNTwst)*bTl) +...
%       mu^3*Axx^2*(0.098*lamNf*bTl^3 + 0.048*thtaN*bTl^4 + 0.049*(thta1+eNTwst)*bTl^5 - 0.023*CmPrfl*cHMR/aLift/cGA*bTl^4);
% 
