function [sys,x0,str,ts] = HorizStabFrcMmtSFunc(t,x,u,flag)
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
% This is the S-function wrapper for the horizontal stabilizer modul
% 
%
%
switch flag,
  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  % Initialize the states, sample times, and state ordering strings.
  case 0
    [sys,x0,str,ts]=mdlInitializeSizes;

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  % Return the outputs of the S-function block.
  case 3
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%
  % Unhandled flags %
  %%%%%%%%%%%%%%%%%%%
  % There are no termination tasks (flag=9) to be handled.
  % Also, there are no continuous or discrete states,
  % so flags 1,2, and 4 are not used, so return an emptyu
  % matrix 
  case { 1, 2, 4, 9 }
    sys=[];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Unexpected flags (error handling)%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Return an error message for unhandled flag values.
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end

% end 
%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts] = mdlInitializeSizes()

sizes = simsizes;
sizes.NumContStates  =  0;
sizes.NumDiscStates  =  0;
sizes.NumOutputs     =  7;  % dynamically sized
sizes.NumInputs      =  9;  % dynamically sized
sizes.DirFeedthrough =  1;  % has direct feedthrough
sizes.NumSampleTimes =  1;

sys = simsizes(sizes);
str = [];
x0  = [];
ts  = [-1 0];   % inherited sample time

% end mdlInitializeSizes

%
%=============================================================================
% mdlOutputs
% Return the output vector of the S-function
%=============================================================================
%
function sys = mdlOutputs(t,x,u)
%----------------------------INPUTS----------------------------------
 vB =[u(1) u(2) u(3)];
 omB=[u(4) u(5) u(6)];

 dle  = u(7);
 lamI = u(8);
 oM   = u(9);
 
 [fTR, mTR, wH] = HorizStabFrcMmtFunc(dle, vB, omB, lamI, oM);
 
 % ----------- OUTPUT --------------
 sys(1) = fTR(1); sys(2) = fTR(2); sys(3) = fTR(3);
 sys(4) = mTR(1); sys(5) = mTR(2); sys(6) = mTR(3);
 sys(7) = wH;
 

% end mdlOutputs

