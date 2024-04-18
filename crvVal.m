function y = crvVal(x,cDat)
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
% This modul calculates values from plynomial coefficients 
% used e.g. to describe fuselage forces and moments. A boundary
% check is performed and linear extrapolation is used outside 
% the boundary values. Uses data read by cXstrct.
%
% Input:
% x: data point at which function value is to be calculated
% cDat: polynomial coefficients
% 
% Output:
% y: function value
%
  if (x > cDat.x2)
    y = cDat.a2*x + cDat.b2;
  elseif (x < cDat.x1)
    y = cDat.a1*x + cDat.b1;
  else
    y = evlPly(x,cDat);
  end
  
  
    y = y*cDat.cAmp;
return
