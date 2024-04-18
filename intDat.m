function out = intDat(y,x,m,idxRe,reFact)
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
% two dimensional linear interpolation in an array y. The interpolation in
% the direction of Reynolds numbers is performed outside and only a ratio
% reFact is handed down into this routine
%
% Input:
% lamNf = velocity component perpendicular to no feathering 
%         (control) plane devided by rotor tip speed (oM*rRot)
% thtaN = collective pitch at blade root
%
% Output:
% value of eq13
% 

   dAlfa= y(m+1,1)-y(m,1);
   xh  = (y(m+1,1)-x)/dAlfa;

   dY  = y(m+1,idxRe) - y(m,idxRe);
   x2  = y(m+1,idxRe) - dY*xh;
   
   idxRe= idxRe-1;
   dY  = y(m+1,idxRe) - y(m,idxRe);
   x1  = y(m+1,idxRe) - dY*xh;
   
   out =x1 + reFact*(x2-x1);