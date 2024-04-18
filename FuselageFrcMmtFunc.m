 function [fFS, mFS] = FuselageFrcMmtFunc( vB, omB, anglC)
%
% Copyright 2010 Juergen Humt
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
%         --- Version 1.0 ---
%
%
% This is the fuselage forces and moments modul.
% The numbers of some of the formulae given in 
% the report are listed in brackets in the code 
% lines
%
% Input:
% vB : velocities in body coordinates
% omB: angular velocities in body coordinates
%
% Output
% fFS : fuselage forces
% mFS : fuselage moments
% 
%
%
 global  rhoAir fpA afH afV lTR hTR ClaFus ClFus CyFus CdFus CmFus...
         detafdalpfh CnFus clNfus cmNfus aDiskMR vfv1

  
  uXB = vB(1);
  uYB = vB(2);
  uZB = vB(3); 
  mu  = anglC(7);
  q = omB(2);
  
  
  % Below influence of downwash on fuselage is calculated. Using clcVix with a factor 
  % of 1.0 gives a highly unstable model. Therefor the wake angle chi is used to scale 
  % vIx. Nevertheless the below formula is pretty much by guess and by golly. The
  % denominator is introduced to make vFx near 1.0 clcVix for mu=0 and slightly larger
  % for mu > 0
  tanChi= atan(uZB/uXB)/(1.57 - 2*mu + mu*mu);
  
  vFx = clcVix(0.25*lTR, 0.5*hTR, vB, anglC);
  xh=rad2deg(atan(vFx/uXB));
  
  vFx=(1 + vfv1*tanChi)*vFx;
  uZB = uZB + vFx; 
  vBXZ2= (uXB^2 + uZB^2);
  vBXY2= (uXB^2 + uYB^2);
  
  % the fuselag anles of attack calculated below are
  % used to calculate coefficient values for cfvVal
  [a36  fhi alfFus]  = alfaFin36(uXB, uZB);
  [b36 fhiB btaFus ] = alfaFin36(uXB, uYB);
  
  % fuselage aerodynamics
  % cd = crvVal(alfFus, CdFus);
  XF = -0.5*rhoAir*vBXZ2*crvVal(alfFus, CdFus)*fpA;  
  YF =  0.5*rhoAir*vBXY2*crvVal(btaFus, CyFus)*afV;
    % rotor download on fuselage added
  ZF = -0.5*rhoAir*clNfus*(vBXZ2*crvVal(alfFus, ClFus) - detafdalpfh*(1-4*mu*mu)*vFx*vFx)*afH ;
  
  % 485.64
  MF=  0.5*rhoAir*cmNfus*((vBXZ2*crvVal(alfFus, CmFus)*fpA)); %  + 1.0*vFx^2*0.8*afH*0.5*lTR) - 1000*afH*(q*0.6667*lTR)^2);
  NF= -0.5*rhoAir*vBXY2*crvVal(btaFus, CnFus)*afV;
  
 
  fFS(1)=XF;
  fFS(2)=YF;
  fFS(3)=ZF;

  mFS(1)= 0;
  mFS(2)= MF;
  mFS(3)= NF;
  
