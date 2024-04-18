 function [fTR, mTR] = TailRotorFrcMmtFunc(dlp, vB, omB)
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
% This is the tail rotor forces and moments modul.
% the numbers of some of the formulae given in the
% report are listed in brackets in the code line
%
% Input:
% dlp: pilot pedal input
% vB : velocities in body coordinates
% omB: angular velocities in body coordinates
% 
% Output
% fTR : tail rotor forces
% mTR : tail rotor moments
%
%
%
 global C6 C7 C8 T1 T2 T3 T4 T5 lTR hTR
  thtaTR=C6*dlp+C7;
  
  if (thtaTR > 0.0873)
     thta2=abs(thtaTR);
  else
  	 thta2= C8;
  end
 
  uXB = vB(1);
  uYB = vB(2);
  
  pB = omB(1);
  rB = omB(3);
 
  vT= uYB - rB*lTR + pB*hTR; % (38) p.7
  
  
  cTR=1.0;
  YTR= cTR*((-T2 + sqrt(T2^2+T3*abs(thtaTR)))^2*thtaTR/abs(1.0e-5+thtaTR)*(1.0+T1*abs(uXB)-...
      (T4+T5*abs(uXB)))/(1.0e-5 + thta2))*vT;  % (39) p.8
 
  % Tail rotor contribution to body forces and moments
  LTR= YTR*hTR;  % (40) p.8
  NTR=-YTR*lTR;  % (41) p.8

  fTR(1)=0;
  fTR(2)=YTR;
  fTR(3)=0;

  mTR(1)= LTR;
  mTR(2)= 0;
  mTR(3)= NTR;
  
  % fTR=[0 0 0];   mTR=[0 0 0];
  