function cT2sa = eq6(mu, lamN, thtaN)
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
% this modul caculates the a rotor parameter of twice the thrust
% coefficient devided by rotor solidity and blade lift coefficient
%
% Input:
% mu    = advance ratio  -> V*cos(alfaD)/oMR
% lamNf = velocity component perpendicular to no feathering 
%         (control) plane devided by rotor tip speed (oM*rRot)
% thtaN = collective pitch at blade root
%
% Output:
% cT2sa  = 0.5*cT/(*sigma*aLift)    see modul description
%

global thta1 bTl gamRot cLftGlo

t31 = 1/2*bTl^2+1/4*mu^2+(1/4/bTl^2-5/36*gamRot^2*bTl^6/(144+gamRot^2*bTl^8))*mu^4;
t32 = 1/3*bTl^3+1/2*bTl*mu^2-4/9/pi*mu^3+(1/3/bTl-25/144*gamRot^2*bTl^7/(144+gamRot^2*bTl^8))*mu^4;
t33 = 1/4*bTl^4+1/4*bTl^2*mu^2+(7/32-2/15*gamRot^2*bTl^8/(144+gamRot^2*bTl^8))*mu^4;

cT2sa = (t31*lamN + t32*thtaN + t33*thta1);