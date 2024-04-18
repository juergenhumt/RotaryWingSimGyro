function out = d2ABCP(dle,dla)
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
% This modul calculates swash plate angles from longitudinal 
% and lateral cyclic stick input
%
% The numbers of some of the formulae given in the report are 
% listed in brackets in the code line
%
% Input:
% dle: longitudinal cyclic stick input
% dla: lateral cyclic stick input
%
% Output
% out = [A1CP, B1CP];
% 
% A1CP: lateral swash plate angle
% B1CP: longitudinal swash plate angle
% 
%
%
  global C1 C4
  fiP=deg2rad(0);
  A1CP=C4*dla*cos(fiP) - C1*dle*sin(fiP); % (23) s.12
  B1CP=C4*dla*sin(fiP) + C1*dle*cos(fiP); % (24) s.12
  % fprintf('A1 %10.4f    B1 %10.4f\n',A1CP,B1CP)
  out = [A1CP, B1CP];
  