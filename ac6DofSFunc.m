%
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
% This is the S function wrapper for the six degree of freedom (dof) body equations
% 
%
%
function [sys,x0,str,ts] = ac6Dof(t,x,u,flag)
% global mHeli, Ixx, Iyy, Izz
% The following outlines the general structure of an S-function.
switch flag

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
    
  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  case { 1, 2, 4, 9 }
      sys=[];
    
  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end
% end sfuntmpl

%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
function [sys,x0,str,ts]=mdlInitializeSizes

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 12;
sizes.NumInputs      = 15;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [-1 0];
% end mdlInitializeSizes

%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
function sys=mdlOutputs(t,x,u)

 [uE, vE, wE, phiDot, thetaDot, psiDot, pDot, qDot, rDot, uDotB, vDotB, wDotB] = ac6DofFunc(u);

%Outputs:
sys(1)  = uE;       % Translatory velocities in EF
sys(2)  = vE;       %
sys(3)  = wE;       %
sys(4)  = phiDot;    % Euler rates
sys(5)  = thetaDot;  %
sys(6)  = psiDot;    %
sys(7)  = pDot;      % Angular accelerations (pDotot, qDotot, rDotot)
sys(8)  = qDot;      %
sys(9)  = rDot;      %
sys(10) = uDotB;     % Translatory acceleration in BF
sys(11) = vDotB;     %
sys(12) = wDotB;     %

% end mdlOutputs
