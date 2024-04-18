function b1 = eq3(mu, lamN, thtaN)
global thta1 wBlade Ib oM gamRot rRot bTl...
       t17 t18 t19 t11N
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
% this modul caculates the rotor lateral 
% flapping angle as per naca 716
%
% Input:
% mu    = advance ratio  -> V*cos(alfaD)/oMR
% lamNf = velocity component perpendicular to no feathering 
%         (control) plane devided by rotor tip speed (oM*rRot)
% thtaN = collective pitch at blade root
%
% Output:
% b1  = rotor lateral flapping angle
%
t17 = 2/9*bTl^2*mu-1/9*mu^3+0.0767/bTl*mu^4+(32/3+7/162*gamRot^2*bTl^8)/(144+gamRot^2*bTl^8)*mu^3;
t18 = 1/6*bTl^3*mu+1/12*bTl*mu^3+0.0175*mu^4+(92/9*bTl+7/216*gamRot^2*bTl^9)/(144+gamRot^2*bTl^9)*mu^3;
t19 = 2/15*bTl^4*mu+2/45*bTl^2*mu^3+0.014*bTl*mu^4+(8*bTl^2+7/270*gamRot^2*bTl^10)/(144+gamRot^2*bTl^10)*mu^3;
t11N =-4/3/bTl*mu+2/3/bTl^3*mu^3-0.14/bTl^4*mu^4;


b1 = gamRot*(t17*lamN+t18*thtaN+t19*thta1); % +t11N*wBlade*0.5*rRot/Ib/oM^2);
