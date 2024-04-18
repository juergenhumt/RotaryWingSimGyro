function cQ2s=  eq9a (mu, lamNf, thtaN, oM)
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
% this modul caculates the accelerating rotor torque
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
% cQ2s  = accellerating rotor torque
global thta1 aLift gamRot bTl cLftGlo cMIN

mu2= mu*mu;
mu4= mu2*mu2;


t41 = 1/2*bTl^2+(5/4+1/1296*gamRot^2*bTl^8)*mu2+...
			(1/2/bTl^2+gamRot^2*bTl^6*(4/3-37/162*gamRot^2*bTl^8-77/46656*gamRot^4*bTl^16)/(144+gamRot^2*bTl^8)^2)*mu4;

t42 = 1/3*bTl^3+(8/3*bTl+1/864*gamRot^2*bTl^9)*mu2+2/9/pi*mu^3+...
			(8/3/bTl+gamRot^2*bTl^7*(224/9-11/648*gamRot^2*bTl^8-41/31104*gamRot^4*bTl^16)/(144+gamRot^2*bTl^8)^2)*mu4;

t43 = 1/4*bTl^4+(2*bTl^2+1/1080*gamRot^2*bTl^10)*mu2+...
			(65/32+gamRot^2*bTl^8*(236/15-7/108*gamRot^2*bTl^8-47/38880*gamRot^4*bTl^16)/(144+gamRot^2*bTl^8)^2)*mu4;

t44 = (8/9*bTl^2+1/2304*gamRot^2*bTl^10)*mu2+...
			(20/9+gamRot^2*bTl^8*(305/36-65/1296*gamRot^2*bTl^8-5/82944*gamRot^4*bTl^16)/(144+gamRot^2*bTl^8)^2)*mu4;

t45 = (4/3*bTl^3+1/1440*gamRot^2*bTl^11)*mu2+...
			(10/3*bTl+gamRot^2*bTl^9*(57/5+7/144*gamRot^2*bTl^8-11/51840*gamRot^4*bTl^16)/(144+gamRot^2*bTl^8)^2)*mu4;

t46 = (1/2*bTl^4+1/3600*gamRot^2*bTl^12)*mu2+...
			(5/4*bTl^2+gamRot^2*bTl^10*(92/25+1/150*gamRot^2*bTl^8-17/129600*gamRot^4*bTl^16)/(144+gamRot^2*bTl^8)^2)*mu4;
			
t47 = -1/108*bTl^5*mu2+(47/9*bTl^8+7/486*gamRot^2*bTl^11)/(144+gamRot^2*bTl^8)*mu4;

t48 = -1/144*bTl^6*mu2+(485/108*bTl^4+5/1296*gamRot^2*bTl^12)/(144+gamRot^2*bTl^8)*mu4;
t49 = -1/180*bTl^7*mu2+(18/5*bTl^5+13/3240*gamRot^2*bTl^13)/(144+gamRot^2*bTl^8)*mu4;

t41N = 1/36*bTl^2*mu2+7/144*mu4;



% cMI = cMIN/oM^2;

cQ2s  =(t41*lamNf^2+t42*lamNf*thtaN+t43*lamNf*thta1+t44*thtaN^2+t45*thtaN*thta1+t46*thta1^2);
% cQ2s  = cQ2s + gamRot*cMI*(t47*lamNf + t48*thtaN + t49*thta1 + t41N*cMI);

cQ2s= aLift*cQ2s;