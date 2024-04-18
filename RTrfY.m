function vOut = RTrfY(thta, uE, vE, wE)
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
% This modul implements a rotation about the body Y axis.
% It is currently used to take into account the main rotor 
% rigging angle. For the inverse transformation call the
% module with an angle of opposite sign
%
% Input:
%   thta: angle about Y axis
%   uE  : x-component of vector to be transformed
%   vE  : y-component of vector to be transformed
%   wE  : z-component of vector to be transformed

% Output:
%  vOut : transformed vector
% 
%
   R_pitch = [...
       cos(thta), 0,-sin(thta);...
            0,    1,      0;...
       sin(thta), 0, cos(thta)];

	  % compute output
  vOut = R_pitch*[uE; vE; wE];
 