function readDmy (inpFile,jRead)
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
% This module performes dummy read actions on a file
% to advance the filepointer b jRead lines, in other
% words jRead lines are skipped in the input file
%
% Input
%  inpFile  input file pointer
%  jRead    number of lines to be skipped
% 
  k=0;
  while ((feof(inpFile) == 0) & k < jRead)
    sAux = fgetl(inpFile);
    k=k+1;
  end 
