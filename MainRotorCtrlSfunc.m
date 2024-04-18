function [sys,x0,str,ts] = MainRotorCtrlSFunc(t,x,u,flag)
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
% This is the S-function wrapper for the main rotor modul
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
    sys=real(mdlOutputs(t,x,u));

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
sizes.NumOutputs     = 18;  % dynamically sized
sizes.NumInputs      = 11;  % dynamically sized
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
% Return the output vector for the S-function
%=============================================================================
%
function sys = mdlOutputs(t,x,u)
%----------------------------INPUTS----------------------------------
 global firstOut nxT jTOut thta1 jTOutEnd
 xDotB = [u(1) u(2) u(3)];
 omDotB= [u(4) u(5) u(6)];
 A1S   = u(7);
 B1S   = u(8);
 thtaN = u(9);
 TMR   = u(10);
 oM    = u(11);
 
 
 xDotC = [0 0 0];
 omDotC= [0 0 0];
 anglC = [0 0 0 0 0];
  
 if abs(TMR) < 1.0e-4
    xh=TMR;
    fprintf('TMR % 17.8e    t %17.8e\n',xh,t);
 end

 [lamNf, xDotC, omDotC, anglC] = MainRotorCntrlFunc(xDotB, omDotB, oM, thtaN, A1S, B1S, TMR);
 %        1        2         3       4        5         6         7         8        9       10       11        12
 sys(1) = xDotC(1);  sys(2) = xDotC(2);  sys(3) = xDotC(3);
 sys(4) = omDotC(1); sys(5) = omDotC(2); sys(6) = omDotC(3);
 sys(7)  =  anglC(1);  sys(8)  =  anglC(2);  sys(9)  =  anglC(3); sys(10) =  anglC(4); 
 sys(11) =  anglC(5);  sys(12) =  anglC(6);  sys(13) =  anglC(7); sys(14) =  anglC(8);
 sys(15) =  anglC(9); sys(16) =  anglC(10);  sys(17) =  anglC(11);
 sys(18) =  lamNf;


 
%  if (t > 1.0e-6)
%    if (firstOut > 0.5)
%      fprintf('a1 %6.2f   A1S %6.2f     b1 %6.2f     B1S %6.2f      thtaN %6.2f\n',rad2deg(anglC(1)),rad2deg(anglC(2)),rad2deg(anglC(3)),rad2deg(anglC(4)),rad2deg(anglC(5)));  
%      firstOut=0;
%      fprintf('jTOut %6.2f  nxT %6.2f\n',jTOut,nxT(jTOut));
%    else
%      if (t > nxT(jTOut)) & (t < nxT(jTOutEnd))
%         fprintf('t %9.6f   a1 %6.2f   A1S %6.2f     b1 %10.7f     B1S %6.2f      thtaN %6.2f      lamNf %15.10f      wC %15.10f\n',t,rad2deg(anglC(1)),rad2deg(anglC(2)),rad2deg(anglC(3)),rad2deg(anglC(4)),rad2deg(anglC(5)),lamNf,xDotC(3));  
%         jTOut=jTOut+1;
%      end
%    end
%  end 

% end mdlOutputs

