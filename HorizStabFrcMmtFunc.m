 function [fHS, mHS, wH] = HorizStabFrcMmtFunc(dle, vB, omB, lamI, oM, anglC)
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
% This is the horizontal stabilizer forces and moments modul.
% the numbers of some of the formulae given in the report are 
% listed in brackets in the code line. In the h-stab coordinate
% system a positive force is up (towards the main rotor) and
% generates a negative pitch moment contribution. A positive
% rate of rotation q thus generates a negative pitch moment.
%
% iRigHS = const -> rigid H-stab
% iRigHS = 1.5e-6-> gearing curve used;
%
% Input:
% dle : longitudinal stick input
% vB  : velocities in body coordinates
% omB : angular velocities in body coordinates
% lamI: induced rotor inflow 
% oM  : rotor rotational velocity rad/s 
%
% Output
% fHS : horizontal stabilizer forces
% mHS : horizontal stabilizer moments
% wH  : test output of velocity resultant
%       in body z direction
% 
%
%
 global  vhv1 H1 H2 H4 hMR lHS hHS areaHS bHS rhoAir iRigHS...
     firstOut firstOutB firstOutC rRot aLiftHS
 
%
% Horizontal stabilizer aerodynamics. Calculate
% H stab angle from longitudinal input.
  R7 = 1/(oM*rRot);

  dls=dle2dlsClc(dle);
  

  vOut = RTrfY(iRigHS+dls, vB(1),  vB(2),  vB(3));
  vIx  = clcVix(lHS, hHS, vB, anglC);
  
  uH = vOut(1);
 % positive vB(3) means, that H-stab moves downward relative to surrounding air
 % induced velocity comes from above so vIx term below has opposite sign of vB(3)
 % positive rotation about y means, that H-stab also moves downward therefore vB(3) 
 % and rotation term -> lHS*omB(2), have the same sign. vIx is multiplied by rotor
 % downwash ratio vhv1 which is zero if HS is outside rotor disk
  wH = vOut(3) + lHS*omB(2) - vhv1*vIx;
  [alHS, phiHS, xh] = alfaFin36(uH,wH);

  % Currently only angles less than 90 degrees are used. Tests on angles 
  % larger than this have to be performed before using the module for this 
  % case
  % For angles larger than 20 degrees the horizontal stabilizer force
  % increases only very slowly if the formulae of naca 73254 are used. 
  % The formula below ensures that force increases nearly linearly 
  % beyond 20 degrees. With this constant slope the q rate is very near 
  % the curve found in naca 73254
  cWH  = 30.00;
  cWH2 =  0.04;
  
  if (abs(xh) < 20)
  % this also includes 160 < alfa < 200
  % uF*vF ca uF^2*tan(alHS)
     ZHS = -H1*wH*uH - cWH*H4*wH*(abs(wH) + cWH2*wH*wH);   % (57) p.9 mod Hu 01/11
%     if firstOut > 0.5
%       fprintf('#1 phiHS %6.2f   xh %6.2f   ZHS %10.2f\n',phiHS,xh,ZHS);  
%       firstOut=0;
%     end 

  %  ZHS = -wH*uH*(H1 + H4*abs(wH)/uH);
  else
    % mod Hu: a factor of 4.15 has been added to smooth
    % out a jump that occured on transition from the
    % 20 degree sector to larger angles
      if ( alHS > 200)
       ZHS=   4.15*H2*uH*uH - cWH*H4*wH*(abs(wH) + cWH2*wH*wH);
    else
       ZHS=  -4.15*H2*uH*uH - cWH*H4*wH*(abs(wH) + cWH2*wH*wH);
    end
    
%     if firstOutB > 0.5
%       fprintf('#2 phiHS %6.2f   xh %6.2f   ZHS %10.2f\n',phiHS,xh,ZHS);  
%       firstOutB=0;
%     end 
  end
  
  % mod Hu:
  % h stab drag is not included in naca 73254
  % aspect ratio of h stab is 4.4 thus lift
  % curve slope is 5.0 (using 5.6 as a0)
  % efficiency factor e=0.9 

  cL = aLiftHS*phiHS;
  cD = 0.1 + cL^2/(pi*bHS*0.9);
  
  fX = 0.5*rhoAir*cD*uH*abs(uH)*areaHS*cos(phiHS)^2/(1+(0.1*wH)^2);

  fHS(1)= -fX;
  fHS(2)=  0;
  fHS(3)= ZHS; % *abs(cos(phiHS));

  mHS(1)= 0;
  mHS(2)= (fHS(3)*lHS - fHS(1)*hHS);
  mHS(3)= 0;

  
  if firstOutC > 0.5
    fprintf('#3 phiHS %6.2f   xh %6.2f   ZHS %10.2f\n',phiHS,xh,ZHS);  
    firstOutC=0;
  end 
