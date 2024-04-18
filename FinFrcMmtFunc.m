 function [fFN, mFN] = FinFrcMmtFunc( vB, omB, dR)
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
% This is the  vertical fin forces and moments modul.
% the numbers of some of the formulae given in the
% report are listed in brackets in the code line
%
% Input:
% vB : velocities in body coordinates
% omB: angular velocities in body coordinates
% 
% Output
% fFN : tail rotor forces
% mFN : tail rotor moments
%
%

 global  rhoAir F1 k1 k2 lVF hVF areaVSfin aLiftVSfin CDovert CDalf9N

 % in forward flight uF > 0
  uF =  vB(1);
  vF = -vB(2) + lVF*omB(3) + 0.5*uF*dR;


  [alFin, phiFin, xh]=alfaFin36(uF,vF);
  
  
  F1 = 0.5*rhoAir*areaVSfin*CDalf9N;  % 0.77
  k1 = 0.5*rhoAir*areaVSfin*aLiftVSfin;    % 1.4

  if (abs(xh) < 40)
  % this also includes 160 < alfa < 200
  % uF*vF ca uF^2*tan(alfaFin)
    YVF= k1*uF*vF + F1*vF*abs(vF);
  else
    k2 = k1*tan(deg2rad(20));      % 0.51
  	if ( alFin > 200)
       YVF= -k2*uF*uF + F1*vF*abs(vF);
    else 
       YVF=  k2*uF*uF + F1*vF*abs(vF);
    end
  end
  
  % fprintf('YVF  %8.3f\n',YVF);
  NVF= -YVF*lVF;  % (53) p.8


  
  fFN(1)= -0.5*rhoAir*CDovert*uF*abs(uF)*areaVSfin*cos(xh);
  fFN(2)= YVF;
  fFN(3)= 0;

  mFN(1)=  fFN(2)*hVF;
  mFN(2)= -fFN(1)*hVF;
  mFN(3)= NVF;
  
 