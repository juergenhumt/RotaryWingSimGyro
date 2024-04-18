function [ret12, untStr] = readLn2(fid)
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
% this moduel reads on line from the input file and converts first token 
% of the input string to a float variable. The rest of the string is 
% scanned for square brackets. If these are found the string enclosed in
% them is returned separately. Usually this is a string describing the 
% unit that has been used for the float value. Default return value (if
% unit string is not present) is "def"
%
% Input
%   fid   input file pointer
%
% Output
%   ret12   float value
%   untStr  unit string, if present
%
  sAux = fgetl(fid);
  [val,rmdr] = strtok(sAux,'%');
  % if there is a closing bracket in the string another unit than standard is given
  idx1=findstr(rmdr,'[')+1;
  if idx1 > 0
    idx2=findstr(rmdr,']')-1;
    if idx2 < 0
      fprintf('error in readLn2, missing delimiter\n');
    else
      untStr =rmdr(idx1:idx2);
    end
  else
     untStr='def';
  end
  
  ret12   = sscanf(val,'%f');
