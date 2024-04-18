function xOut = muCtIt(x, vFlight, thtaN, thtaClmb, cOut)
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
% This modul is used in an iterative search of rotor speed and
% disk angle for a given value of advance ratio and disk angle
%
% Input:
%   x       : vector of input values
%             x(1) tip speed
%             x(2) disk angle
%   vFlight : flight speed
%   thtaN   : collective pitch at blade root
%   thtaClmb: angle of climb/descent
%   cOut    : ouput control parameter
%
% Output:
%   xOut    : vector of output values
%             xOut(1) difference of tip speed required minus tip speed in
%             xOut(2) differenct of disk anlge required minus disk angle in
% 

global rhoAir aDiskMR rRot mHeliG sigma iRigMR epsGlo


oMR   = x(1);
alfaD = x(2);


thtaTrim=alfaD + thtaClmb - iRigMR;


xCTNN = rhoAir*aDiskMR*(oMR)^2;
cTs = mHeliG*cos(thtaTrim)/xCTNN/sigma;

xL=1.0e-4; xR= 0.75; itMax=150;
[muIt, fMin, cTsIt, cTIt, cL, alfaNf, a1MR, oMRIt, vFIt, nRot, nRotMin, lAmNf, k] = gSearchOm(xL,xR,epsGlo,itMax,thtaN, cTs, thtaTrim);
alfaDIt = alfaNf + a1MR;

vF = norm(vFlight);
mu = vF/oMR*cos(alfaD);

if cOut < 0.5
  xOut(1)=mu - muIt;
  xOut(2)=alfaD -alfaDIt;
else
  xOut(1)= muIt;
  xOut(2)= alfaDIt;
  xOut(3)= alfaNf;
  xOut(4)= a1MR;
  xOut(5)= thtaTrim;
  xOut(6)= cTs;
end