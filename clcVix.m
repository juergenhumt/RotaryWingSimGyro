function vIx = clcVix(lX, hX, vB, anglC)
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
% induced velocity (rotor downwash) at the location lX hX  is calculated
% using the values published by Heyson and Katzoff. This routine is used 
% in calculation of effective angle of attack for components in rotor 
% downwash: horizontal, stabilizer, wing and fuselage
%
% Input:
% lX: logitudinal distance from rotor hub to point where downwash is calculated
% hX: vertical distance from rotor hub to point where downwash is calculated
% vB: velocity in body coordinates
% anglC: rotor angles
%
% Output:
% vIx: induced velocity
%
global xCG hMR iRigMR rRot


  lamI  = anglC(6);
  % zero coordinate of model is at CoG but induced velocity
  % calculation is based on distance from rotor hub 
  ksi = lX/rRot;  % ksi   = (lX-xCG)/rRot;
  oMR   = anglC(8)*rRot;
  
  % the distance between hub plane and HS is:
  zeta = hX/rRot; % zeta  = (hX - hMR)/rRot;
  
  % epsN is the angle between free stream velocity 
  % and the vortex sheets emanating from the rotor
  alfaD = anglC(12) + anglC(1);
  epsN  = lamI/(sqrt(vB(1)^2 + vB(3)^2)*cos(alfaD)/oMR) ;
  
  % alfHT is the angle between hubplane and the
  % vortex sheets emanating from the rotor
  alfHT = anglC(13) + anglC(2) - iRigMR - epsN;
  
  % the downward distance the sheets have traveled
  % when reaching the HS is in nondimensional form
  zetaS= ksi*alfHT;
  
  % the distance between the point where the vortex 
  % sheet springs from the rotor disk and HS is
  zetaR = zetaS - zeta;
  
  % the amplification factor depending on the relative
  % position of HS in the flow field created by the 
  % induced flow is:
  cZeta = 0.8*zetaIntp(ksi,zetaR);
  
  % the induced velocity in the flow field at that point is
  vIx = lamI*oMR*cZeta;

