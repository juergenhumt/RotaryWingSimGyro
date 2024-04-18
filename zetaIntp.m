function z=zetaIntp(ksi,zeta)
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
% This module interpolates between the two curves for induced velocity 
% amplification factor given by Bramwell
%
% Input
%   ksi  nondimensional longitudinal distance to rotor hub (distance devided by rotor radius)
%   zeta nondimensional vertical distance to rotor hub (distance devided by rotor radius)
% 
% Output 
%   z  amplification factor
% 
global zPA zPB
zA=polyval(zPA,zeta);
zB=polyval(zPB,zeta);

% ksi = 1.07-> z=zA
% ksi = 2.07-> z=zB

z=(2.07-ksi)/(2.07-1.07)*(zA-zB)+zB;

