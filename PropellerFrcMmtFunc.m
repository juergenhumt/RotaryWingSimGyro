 function [fPR, mPR] = PropellerFrcMmtFunc( vB, thrtlStng)
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
% This is the propeller forces and moments modul.
%
% Input:
%   vB       : velocities in body coordinates
%   thrtlStng: throttle setting >=0 (if greater one (=1.0)
%              gyro uses more power than installed)
%
% Output
%   fPR : propeller forces
%   mPR : propeller moments
%
% Global:
%   hPR : vertical distance of propeller line of action
%         from center of gravity( i.e. hPR was rotated to
%         to be perpendicular to the prop thrust line)
%   lPR : analogus 
% tauProp : angle at which propeller is inclined. If greater
%           zero propeller points downward for a tractor
%
 global  hPR lPR tauProp

% body coordinates are rotated to the prop system
  vOut= RTrfY(tauProp,vB(1),vB(2),vB(3));

% propeller thrust
  XPR= thrustCalcUS(vOut(1), thrtlStng);
%   if XPR < 1.0
%     XPR=XPR*XPR;
%   end
% rotate forces and moments to body coordinates
  fPR = RTrfY(-tauProp,XPR, 0,0)';
%'  
  mY  = -(fPR(1)*hPR - fPR(3)*lPR);
  mX  = 0.0; % 7.0e-2*mY;
  mPR = [mX mY 0];

