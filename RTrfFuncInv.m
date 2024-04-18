function vOut = RTrfFuncInv(phi, thtaInp, psi, ue, ve, we)
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
% This modul implements the  body to local (earth) frame 
% transformation matrices found in naca cr 2497. Since
% the body frame is along the shaft, the main rotor
% rigging angle has to be added to the pitch component
% 
% Input:
%   phi : roll angle
%   thta: pitch angle
%   psi : yaw angle
%   uE  : x-component of vector to be transformed
%   vE  : y-component of vector to be transformed
%   wE  : z-component of vector to be transformed

% Output:
%  vOut : transformed vector
% 
%
% 
global iRigMR
      R_roll = [...
	          1,     0,        0;...
	          0,  cos(phi), sin(phi);...
	          0, -sin(phi), cos(phi)];
      thta=thtaInp; % + iRigMR;
	  R_pitch = [...
	          cos(thta), 0,-sin(thta);...
	               0,    1,      0;...
	          sin(thta), 0, cos(thta)];

	  R_yaw = [...
	           cos(psi), sin(psi), 0;...
	          -sin(psi), cos(psi), 0;...
	              0,         0,    1];

	  % compute output
	  Rb2e = R_yaw'*R_pitch'*R_roll'; %'
	  
	  vOut = Rb2e*[ue; ve; we];
 