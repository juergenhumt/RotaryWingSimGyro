function [alfa, phi, xh]=alfaFin36(uF, vF)
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
%         --- Version 2.0 ---
%
%
% This module calculates the angle of attack for a 
% horizontal plane or vertical fin given two body 
% frame velocities 
% 
% Input:
% uF : velocity along body x direction
% vF : velocity along body y direction
%
% Output
% alfa : angle of attack from x axis having a 
%        value between 0 and 360 degrees
% phi  : angle of attack from x axis having a 
%        value betwee -pi < 0 < +pi (radians)
% xh   : same as phi but in degree       
%
%
% 
  phi=asin(vF/sqrt(uF^2+vF^2));
  xh =rad2deg(phi);

  if (uF < 0);
     alfa = 180 - xh;
  else
    if (vF < 0)
      alfa = 360 + xh;
    else
      alfa = xh;
    end
  end
