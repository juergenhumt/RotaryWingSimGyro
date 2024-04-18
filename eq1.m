function aN = eq1(mu, lamN, thtaN)
global thta1 wBlade Ib oM gamRot rRot bTl...
       t11 t12  t13
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
% this modul caculates the rotor coning angle
% as per naca 716
%
% Input:
% mu    = advance ratio  -> V*cos(alfaD)/oMR
% lamNf = velocity component perpendicular to no feathering 
%         (control) plane devided by rotor tip speed (oM*rRot)
% thtaN = collective pitch at blade root
%
% Output:
% aN  = coning angle
t11 = 1/6*bTl^3+0.040*mu^3-5/144*(gamRot^2*bTl^7/(144+gamRot^2*bTl^8))*mu^4;
t12 = 1/8*bTl^4+1/8*bTl^2*mu^2-1/64*mu^4-25/576*(gamRot^2*bTl^8/(144+gamRot^2*bTl^8))*mu^4;
t13 = 1/10*bTl^5+1/12*bTl^3*mu^2-1/30*(gamRot^2*bTl^9/(144+gamRot^2*bTl^8))*mu^4;

aN = gamRot*(t11*lamN+t12*thtaN+t13*thta1); % -wBlade/Ib/oM^2/gamRot);
