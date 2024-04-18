function cQ2sd =  eq11d (mu, lamN, thtaN, oM)
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
% this modul caculates the decelerating rotor torque
% as per naca 716
%
% Input:
% mu    = advance ratio  -> V*cos(alfaD)/oMR
% lamNf = velocity component perpendicular to no feathering 
%         (control) plane devided by rotor tip speed (oM*rRot)
% thtaN = collective pitch at blade root
% oM    = rotor speed in rad/sec
%
% Output:
% cQ2sd  = decelerating rotor torque
global thta1 sigma gamRot bTl dDrgN dDrg1 dDrg2 cMIN

t51 = 1/4+1/4*mu^2-1/32*mu^4;

t52 = 1/3 -5/72*gamRot^2*bTl^5/(144+gamRot^2*bTl^8)*mu^4;

t53 = 1/4+1/4*mu^2 - 25/288*gamRot^2*bTl^6/(144+gamRot^2*bTl^8)*mu^4;

t54 = 1/5+1/6*mu^2 -1/15*gamRot^2*bTl^7/(144+gamRot^2*bTl^8)*mu^4;

t55 = 1/2+(-1/4+1/bTl^2+1/2/bTl^4+gamRot^2*(1/162*bTl^4-1/81*bTl^5+1/144*bTl^6))*mu^2+...
			(-3/4/bTl^2+1/bTl^4+1/4/bTl^6+gamRot^2*(-1/162*bTl^2+1/162*bTl^3+1/324*bTl^4-1/576*bTl^6)+...
      gamRot^2*(7/9*bTl^2-37/27*bTl^3-13/27*bTl^4)/(144+gamRot^2*bTl^8)+...
      gamRot^2*128*bTl^2/(144+gamRot^2*bTl^8)^2+...
      gamRot^4*(7/2916*bTl^10-7/1458*bTl^11-7/2592*bTl^12)/(144+gamRot^2*bTl^8)^2+...
			gamRot^4*193/162*bTl^10/(144+gamRot^2*bTl^8)^2+...      
			gamRot^6*49/23328*bTl^18/(144+gamRot^2*bTl^8)^2)*mu^4;

t56 = 2/3+(4/3/bTl+4/3/bTl^3+gamRot^2*(1/108*bTl^5-1/54*bTl^6+1/96*bTl^7))*mu^2+4/9/pi*mu^3+...
			(-1/bTl+8/3/bTl^3+1/bTl^5+gamRot^2*(-1/108*bTl^4+13/864*bTl^5-1/384*bTl^7)+...
      gamRot^2*(161/108*bTl^4-811/324*bTl^5-113/108*bTl^7)/(144+gamRot^2*bTl^8)+...
      gamRot^2*736/3*bTl^3/(144+gamRot^2*bTl^8)^2+...
      gamRot^4*(7/1944*bTl^11-7/972*bTl^12-7/1728*bTl^13)/(144+gamRot^2*bTl^8)+...
			gamRot^4*233/108*bTl^11/(144+gamRot^2*bTl^8)^2+...      
			gamRot^6*49/15552*bTl^19/(144+gamRot^2*bTl^8)^2)*mu^4;

t57 = 1/2+(1+1/bTl^2+gamRot^2*(1/135*bTl^6-2/135*bTl^7+1/120*bTl^8))*mu^2+...
			(-11/16+2/bTl^2+3/4/bTl^4+gamRot^2*(-1/810*bTl^4-2/405*bTl^5+23/2160*bTl^6-1/408*bTl^8)+...
      gamRot^2*(157/135*bTl^4-37/18*bTl^5-13/18*bTl^6)/(144+gamRot^2*bTl^8)+...
      gamRot^2*192*bTl^4/(144+gamRot^2*bTl^8)^2+...
      gamRot^4*(7/2430*bTl^12-7/1215*bTl^13-7/2160*bTl^14)/(144+gamRot^2*bTl^8)+...
			gamRot^4*229/135*bTl^12/(144+gamRot^2*bTl^8)^2+...
			gamRot^6*49/19440*bTl^20/(144+gamRot^2*bTl^8)^2)*mu^4;

t58 = 1/4+(1/4+8/9/bTl^2+gamRot^2*(1/288*bTl^6-1/144*bTl^7+1/256*bTl^8))*mu^2+...
			(-1/32+4/3/bTl^2+8/9/bTl^4+gamRot^2*(1/288*bTl^4-1/96*bTl^5+11/1152*bTl^6-1/1024*bTl^8)+...
      gamRot^2*(119/162*bTl^4-94/81*bTl^5-47/72*bTl^6)/(144+gamRot^2*bTl^8)+...
      gamRot^2*1058/9*bTl^4/(144+gamRot^2*bTl^8)^2+...
      gamRot^4*(7/5184*bTl^12-7/2592*bTl^13-7/4608*bTl^14)/(144+gamRot^2*bTl^8)+...
			gamRot^4*2557/2892*bTl^12/(144+gamRot^2*bTl^8)^2+...
			gamRot^6*49/41472*bTl^20/(144+gamRot^2*bTl^8)^2)*mu^4;

t59 = 2/5+(1/3+4/3/bTl+gamRot^2*(1/180*bTl^7-1/90*bTl^8+1/160*bTl^9))*mu^2+...
			(2/bTl+4/3/bTl^3+gamRot^2*(1/216*bTl^5-2/135*bTl^6+41/2880*bTl^6-1/640*bTl^9)+...
      gamRot^2*(617/540*bTl^5-2087/1080*bTl^6-107/120*bTl^7)/(144+gamRot^2*bTl^8)+...
      gamRot^2*184*bTl^5/(144+gamRot^2*bTl^8)^2+...
      gamRot^4*(7/3240*bTl^13-7/1620*bTl^14-7/2880*bTl^15)/(144+gamRot^2*bTl^8)+...
			gamRot^4*31/20*bTl^13/(144+gamRot^2*bTl^8)^2+...
			gamRot^6*49/25920*bTl^21/(144+gamRot^2*bTl^8)^2)*mu^4;

t51N = 1/6+(5/8+gamRot^2*(1/450*bTl^8-1/225*bTl^9+1/400*bTl^10))*mu^2+...
			(3/4+1/2/bTl^2+gamRot^2*(1/675*bTl^6-7/1350*bTl^7+19/3600*bTl^8-1/1600*bTl^10)+...
      gamRot^2*(4/9*bTl^6-4/5*bTl^7-3/10*bTl^7)/(144+gamRot^2*bTl^8)+...
      gamRot^2*72*bTl^6/(144+gamRot^2*bTl^8)^2+...
      gamRot^4*(7/8100*bTl^14-7/4050*bTl^15-7/5400*bTl^16)/(144+gamRot^2*bTl^8)+...
			gamRot^4*137/225*bTl^14/(144+gamRot^2*bTl^8)^2+...
			gamRot^6*49/64800*bTl^21/(144+gamRot^2*bTl^8)^2)*mu^4;

t511 = (-2/27*bTl+4/27*bTl^2-1/12*bTl^3)*mu^2+...
			(-2/27+2/27/bTl-1/27*bTl+1/48*bTl^3+...
      (64/9+32/9/bTl-4*bTl)/(144+gamRot^2*bTl^8)+...
      gamRot^2*(-7/468*bTl^7-7/243*bTl^8+7/432*bTl^9)/(144+gamRot^2*bTl^8))*mu^4;

t512 = (-1/18*bTl^2+1/9*bTl^3-1/16*bTl^4)*mu^2+...
			(1/18*bTl-13/144*bTl^2-1/64*bTl^4+...
      (-92/27+184/27*bTl+23/6*bTl^2)/(144+gamRot^2*bTl^8)+...
      gamRot^2*(-7/648*bTl^8+7/323*bTl^9+7/576*bTl^10)/(144+gamRot^2*bTl^8))*mu^4;

t513 = (-2/45*bTl^3+4/45*bTl^4-1/20*bTl^5)*mu^2+...
			(1/135*bTl+4/135*bTl^2-23/360*bTl^3+1/80*bTl^5+...
      (8/3*bTl+16/3*bTl^2+3*bTl^3)/(144+gamRot^2*bTl^8)+...
      gamRot^2*(-7/810*bTl^9+7/405*bTl^10+7/720*bTl^11)/(144+gamRot^2*bTl^8))*mu^4;

t514 = (1/4-4/9/bTl+2/9/bTl^2)*mu^2+...
       (-1/16+1/9*bTl^2+2/9/bTl^3-2/9/bTl^4)*mu^4;

% coeff716(mu);        
cQ2sd = (dDrgN*t51+dDrg1*(t52*lamN+t53*thtaN+t54*thta1)...
    + dDrg2*(t55*lamN^2+t56*lamN*thtaN+t57*lamN*thta1...
    + t58*thtaN^2+t59*thtaN*thta1+t51N*thta1^2));


% cMI = cMIN/oM^2;
% cQ2sd = cQ2sd + gamRot*cMI*(t511*lamN + 512*thtaN + t513*thta1 + t514*cMI);
