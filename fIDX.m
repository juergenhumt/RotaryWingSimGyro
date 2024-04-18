function jOut=fIDX(y,x)
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
% This module (fIDX = find index) returns the index of an element in a 
% vector y with monotonically increasing first column. It is the index 
% of the element that has the value closest to but below the value of x
%
% Input:
% y: vector with element values increasing montonically 
%
% Output:
% jOut: index 
% 
  j=1; [m,n]=size(y);
  
  k=bitshift(m,-1);
  
  while (k~=j) 
       if (x < y(k))
           m=k;
       else
           j=k;
       end
      k= bitshift(j+m,-1);     
  end 
 jOut=k;