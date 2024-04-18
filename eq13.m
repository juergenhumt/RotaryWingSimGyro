function out =  eq13 (mu,lamNf, thtaN)
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
% This modul calculates equation 13 page 214 as per naca 716
% using the variable names of this program the left hand side
% of equation 13 would be written as:
% mu*2*cT/(sigma*aLift)*DvLN
%
% Input:
% lamNf = velocity component perpendicular to no feathering 
%         (control) plane devided by rotor tip speed (oM*rRot)
% thtaN = collective pitch at blade root
%
% Output:
% out   = value of eq13
% 
global thta1 aLift sigma cMI dDrgN dDrg1 dDrg2...
        t61 t62 t63 t64 t65 t66 t67 t68 t69 t61N...
        t611 t612 t613 t614

        
coeff716(mu);
out = (dDrgN/aLift*t61+dDrg1/aLift*(t62*lamNf+t63*thtaN+t64*thta1)+...
     dDrg2/aLift*(t65*lamNf^2+t66*lamNf*thtaN+t67*lamNf*thta1+t68*thtaN^2+...
     t69*thtaN*thta1+t61N*thta1^2));
