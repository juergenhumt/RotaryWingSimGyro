function [aNr, a1r, a2r, b1r, b2r, dQ] = clcRotAngl(mu, lamNf, thtaN, oM, rotThrst)
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
% This module calculates all rotor angles
% for a give flight state 
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

global rRot

% coning angle aN 
  aNr = eq1(mu, lamNf, thtaN);

% a1 and b1 are rotor flap angles refered to the no feathering axis
  a1r = eq2(mu, lamNf, thtaN);  % angular terms from eq (6)
  b1r = eq3(mu, lamNf, thtaN);  % angular terms from eq (7)
  
  a2r = 0;
  b2r = 0;
  
  
  dQ = eq9a(mu,lamNf,thtaN,oM) - eq11d(mu,lamNf,thtaN,oM);
  
