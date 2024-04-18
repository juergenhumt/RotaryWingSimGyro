function [uE, vE, wE, phiDot, thetaDot, psiDot, pDot, qDot, rDot, uDotB, vDotB, wDotB] = ac6DofFunc(u)
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
% This  modul calculates translatory and rotatory accellerations.
% The formulae in this module are not part of NACA tm 73254. These
% are the usual equations of motion for a 6 degree of freedom 
% aircraft model:
%
% Input:
% see below in code
% 
% Output
% see below in code
%
%

  global mHeli Ixx Iyy Izz
	% Input:
	phi =  u(1);    % Euler angles
	theta= u(2);    %
	psi =  u(3);    %
	fx  =  u(4);    % Force components
	fy  =  u(5);    %
	fz  =  u(6);    %
	p   =  u(7);    % Angular velocities
	q   =  u(8);    %
	r   =  u(9);    %
	L   =  u(10);   % Torques
	M   =  u(11);   %
	N   =  u(12);   %
	ux  =  u(13);   % Time delayed body velocities
	uy  =  u(14);   %
	uz  =  u(15);   %

  % Output:    
  % From the input variables linear and rotational velocites 
  % in earth cooridnates are calculated as well as ....

  % V - translatory velocity transformation to earth frame
    vOut  = RTrfFuncInv(phi, theta, psi, ux, uy, uz);
  % V - translatory velocity components, earth frame
    uE=vOut(1); vE=vOut(2); wE=vOut(3);

  % 
  % ... linear and rotational velocites 
  % in body fixed coordinates

	%ThetaDot - Euler rates
	phiDot   = p+sin(phi)*tan(theta)*q+cos(phi)*tan(theta)*r;
	thetaDot = cos(phi)*q-sin(phi)*r;
	psiDot   = (sin(phi)/cos(theta))*q+(cos(phi)/cos(theta))*r;

	%OmegaDot - Angular accelerations (pDot, qDot, rDot (body frame))
	pDot =  ((Iyy-Izz)*q*r+L)/Ixx;
	qDot =  ((Izz-Ixx)*p*r+M)/Iyy;
	rDot =  ((Ixx-Iyy)*p*q+N)/Izz;

    % pDot=0; qDot=0; rDot=0;
    
	%vDotBot - Translatory acceleration in body frame
	uDotB = 1/mHeli*fx-q*uz+r*uy;
	vDotB = 1/mHeli*fy+p*uz-r*ux;
	wDotB = 1/mHeli*fz-p*uy+q*ux;
    
 