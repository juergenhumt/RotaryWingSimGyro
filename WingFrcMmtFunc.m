 function [fWN, mWN] = WingFrcMmtFunc(vB, omB, oM, anglC)
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
% This is the wing forces and moments modul.
% In the wing coordinate system a positive force is up (towards 
% the main rotor) 
%
% Input:
% vB : velocities in body coordinates
% omB: angular velocities in body coordinates
% lamNf: iterated rotor inflow 
% oM   : rotational speed of rotor (2*pi*n) 
%
% Output
% fWN : wing forces
% mWN : wing moments
% 
%
%
 global rhoAir kinVisc firstOutW logFile...
   lWN hWN hMR xCG aWing bWing cRwing cTwing alpnWing...
   iRigWing rRot aLiftWing CLwngGlo CDwngGlo
 
%
% Wing aerodynamics. Calculate wing angle of attack.

  vqAux = lWN*omB(2);
  vOut = RTrfY(iRigWing + alpnWing, vB(1),  vB(2),  vB(3));
  
  x1= (lWN)/rRot;  x2= (hMR - hWN)/rRot;
  vIx  = clcVix(x1, x2, vOut, anglC);
  
  uH = vOut(1);
  wH = vOut(3) + vqAux; %  + vIx;
  
  [alWN, phiWN, xh]=alfaFin36(uH,wH);

  reLoc= norm(vB)*cRwing/kinVisc;
  
  [cL cDN cM err] = afDataInt(phiWN,reLoc);
  % cL = clcCl(rad2deg(phiWN));
  if cL > 0
    cL = CLwngGlo*cL;
    ZWN =  0.5*rhoAir*cL*aWing*vOut(1)*vOut(1);
  else
    ZWN =  0.5*rhoAir*(uH*uH + wH*wH)*aWing*sin(phiWN)/(1+(0.1*wH)^2);
  end
  
  cDN= CDwngGlo*cDN;
  cD = cDN + cL^2/(pi*bWing*0.9);
  fX =  0.5*rhoAir*cD*uH*abs(uH)*aWing*cos(phiWN)^2/(1+(0.1*wH)^2);
  
  if (abs(xh) < 20)
  % this also includes 160 < alfa < 200

     
    if firstOutW > 0.5
      fprintf('#1 phiWN %6.2f   xh %6.2f   ZWN %10.2f\n',phiWN,reLoc,ZWN);  
      firstOutW=0;
    end 

  else
%     % airflow from trailing to leading edge
%     if ( alWN > 90)
%        fX=-fX;
%     end
%     
%     if ( alWN > 200)
%        ZWN=  ZWN;
%     else
%        ZWN=  -ZWN;
%     end
%     % fprintf('phiWN %6.2f   xh %6.2f   ZWN %10.2f\n',phiWN,xh,ZWN);  
%     
    if firstOutW > 0.5
      fprintf('#2 phiWN %6.2f   xh %6.2f   ZWN %10.2f\n',phiWN,xh,ZWN);  
      firstOutW=0;
    end 
  end

  fwnLoc(1)= -fX;
  fwnLoc(2)=  0;
  fwnLoc(3)= -ZWN;

  fWN = RTrfY(-(iRigWing + alpnWing),fwnLoc(1),fwnLoc(2),fwnLoc(3))';

  % fWN(3)=0;
  
  mWN(1)= 0;
  mWN(2)= (fWN(3)*lWN - fWN(1)*hWN) - 0.5*rhoAir*cM*(cTwing+cRwing)*vOut(1)*vOut(1)*aWing;
  mWN(3)= 0;

  % fprintf(logFile,'pfhi %6.2f cL   %10.5f  fWN(1)  %15.2f   fWN(2)  %15.2f   fWN(3)  %15.2f\n',xh,cL,fWN(1),fWN(2),fWN(3));
   
  if firstOutW > 0.5
     fprintf('#3 phiWN %6.2f   xh %10.5e   ZWN %10.2f\n',phiWN,reLoc,ZWN);  
     firstOutW=0;
  end 
