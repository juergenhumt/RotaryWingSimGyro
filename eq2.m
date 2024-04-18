function a1 = eq2(mu, lamNf, thtaN)
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
% this modul caculates the rotor longitudinal 
% flapping angle a1 as per naca 716
%
% Input:
% mu    = advance ratio  -> V*cos(alfaNf)/oMR
% lamNf = velocity component perpendicular to no feathering 
%         (control) plane (=V*sin(alfaNf)) devided by rotor tip speed (oM*rRot)
% thtaN = collective pitch at blade root
%
% Output:
% a1  = rotor longitudinal flapping angle

global thta1 t14 t15 t16 cMI gamRot bTl
t14 = 2/bTl^2*mu+1/2/bTl^4*mu^3+10/27*(gamRot^2*bTl^4/(144+gamRot^2*bTl^8))*mu^3;
t15 = 8/3/bTl*mu+4/3/bTl^3*mu^3+0.212/bTl^4*mu^4+25/54*(gamRot^2*bTl^5/(144+gamRot^2*bTl^8))*mu^3;
t16 = 2*mu+1/bTl^2*mu^3+16/45*(gamRot^2*bTl^6/(144+gamRot^2*bTl^8))*mu^3;
%                           
a1 = (t14*lamNf + t15*thtaN + t16*thta1);
