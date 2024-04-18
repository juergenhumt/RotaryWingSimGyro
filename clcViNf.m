function [cTs, cT, lamNf, lamIj, kG, jItr]= clcViNf(mu, wNf, oM, thtaN, pC, zB, cTReq)
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
% This modul performs the inflow calculation using Padfields
% iterative scheme
%
% Input:
%   mu    = advance ratio  -> V*cos(alfaD)/oMR
%   wNf    = velocity component perpendicular to no feathering (control) plane 
%   thtaN = collective pitch at blade root
%   pC    = control (swash plate) roll rate
%   zB    = height above ground
%   cTReq = required thrust coefficient
%  
% Output:
%   cTs   = iterated thrust coefficient devided by solidity (sigma)
%   cT    = iterated thrust coefficient
%   lamNf = inflow at no feathering axis (control/swash plate)
%   lamIj = induced inflow
%   kG    = ground effect factor
%   jItr  = number of iterations (usually 5)
%
%
%
global sigma aLift thta1 bTl rRot cLftGlo cTCx clcCT_hndl...

   
stopIt = 100;
epsItr =1.0e-12;
cDmp   =1.0;

hJ=1;
oMR= oM*rRot;
R2=1/(4*oM);

G1  = 0.2;
kG  = 1-exp(-zB/(2*rRot)/G1); % ground effect factor

% induced velocity for hover taken as start value
% lamStart = sqrt(cT/2) = sqrt(tc*sigma/2)
lamIj= sqrt(0.5*cTReq);

b5 = bTl^5;
b2 = bTl^2;

% For inflow calculation the air moving upward through
% the rotor was taken positive and induced velocity was
% negative. This is a left over from the analyis of the
% autogyro. Thus the minus below.
muZ=  -wNf/oMR;

jItr    = 0;
while ((abs(hJ) > epsItr) & (jItr < stopIt))
      
   lamNf = (muZ - lamIj);
   LAMDA=mu^2+(lamNf)^2;
   % T from naca - 73254 the formula in brackets gives 2*cT/(sigma*aLift)
   cT = cTCx*feval(clcCT_hndl,mu,lamNf,thtaN,oM,cTReq)*kG;
   % cT = cT + c1*R2*mu*pC;
   % fprintf('%s%10.5f\n','cT Itr     ',cTReq);
   gN = lamIj - cT/(2*LAMDA^0.5);
   hJ = (2*lamIj*LAMDA^0.5 - cT)*LAMDA/(2*LAMDA^1.5+ 0.5*cTCx*LAMDA - cT*(lamNf));
   lamIj = lamIj - cDmp*hJ;
   jItr=jItr+1;
end

if jItr==stopIt
    fprintf('Warning jItr=%d\n',jItr);
end

cTs=cT/sigma;

