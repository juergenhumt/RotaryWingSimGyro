function out=dle2dlsClc(dle)
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
% This module calculates the horizontal stabilizer angle 
% for a given longitudinal cyclic input by evaluating
% the polynomial for the gearing curve
% 
% Input:
% dle : longitudinal cyclic stick motion
% 
% Output
% horizontal stabilizer angle
%
%
% 
global dle2dlsP dleOffs dleMax dleMin fiHsMin fiHsMax iRigMR

   if dle < dleMin
     out=fiHsMin;
   else 
     if dle > dleMax
       out=fiHsMax;
     else
	   out = polyval(dle2dlsP,dle);
	 end
   end
   
   out = out  + dleOffs;
