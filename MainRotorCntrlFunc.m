function [lamNfOut, vCt, omCt, anglC] = MainRotorCntrlFunc(vBd, omBd, oM, thtaN, A1S, B1S, TMR)
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
% The numbers of some of the formulae given in the
% report are listed in brackets in the code line
%
% This function calculates rotor inflow and control 
% angles from translational and rotational velocities
% and pilot control inputs.
%
% Input:
% vB     : velocities in body coordinates
% omB    : angular velocities in body coordinates
% oM     : rotor speed in rad/s
% thtaN  : pilot collective pitch input
% A1S,B1S: pilot longitudinal and lateral swash plate input
%
% Output:
% lamNfOut: inflow calculated using Padfields iterative scheme
% vCt   : velocities in control plane (swash plate) coordinates
% omCt  : angular velocities in control plane (swash plate) coordinates
% anglC: control angles
% anglC(1) = a1   : first harmonic longitudinal flapping coefficient
% anglC(2) = A1S  : pilot lateral swash plate input
% anglC(3) = b1   : first harmonic  lateral flapping coefficient
% anglC(4) = B1S  : pilot longitudinal swash plate input
% anglC(5) = thtaN: collective pitch angle 
% 
%
%
  global...
    rhoAir aN kB hMR thta1 aDiskMR mHeliG rRot iRigMR gamRot clcRtAngl_hndl


  oMR = oM*rRot;
  % fprintf('\nnaca 73254 control\n')
  xCTNN=rhoAir*aDiskMR*(oMR)^2;
  cT= TMR/xCTNN;
  
  s=0;
  
%   uB= vBd(1);  % longitudinal velocity, body axis
%   vB= vBd(2);  % lateral velocity, body axis
%   wB= vBd(3);  % vertical velocity, body axis
% 
%   pB= omBd(1);  % body axis roll rate
%   qB= omBd(2);  % body axis pitch rate
%   rB= omBd(3);
%   
  vOut= RTrfY( iRigMR,vBd(1),vBd(2),vBd(3));

  uB= vOut(1); vB= vOut(2); wB= vOut(3);

  vOut= RTrfY( iRigMR,omBd(1),omBd(2),omBd(3));
  pB= vOut(1); qB= vOut(2); rB= vOut(3);
  %   
  % amplitude and azimuth of maximum swash plate
  % deflection from combined action of longitudinal
  % and lateral swash plate Input
  % Hu replaced A/B1C by A/B1S
  % fprintf('A1S  %10.5f        B1S %10.5f\n',A1S,B1S);
  fiC=acos(cos(A1S)*cos(B1S));
  xh=sin(B1S);
  if abs(xh) > 1.0e-10
    fiT=atan(tan(A1S)/xh);
  % modified calculation of longitudinal and lateral cyclic
  % in hub wind axes. Added negative branch of if statement 
  % for small B1S.
  elseif xh < 0
    fiT= -0.5*pi;
  else
    fiT=  0.5*pi;
  end 	
 	

  % uC: longitudinal velocity in wind control axis
  uC= -uB*cos(B1S)-wB*sin(B1S) + qB*hMR;  % (11) s.11
  % vC: lateral velocity in wind control axis
  vC= -vB*cos(A1S)-wB*sin(A1S) - pB*hMR;  % (12) s.11
  % wC: vertical velocity, wind control axis
  wC= -wB*cos(A1S)*cos(B1S) + sin(B1S)*(uB-qB*hMR) + sin(A1S)*(vB+ pB*hMR); % (13) s.11
 
  xNN = sqrt(uC^2+vC^2);
  mu = xNN/oMR;
  % coeff716(mu);
  tauR = 16/(gamRot*oM);  % 0.144 main rotor time constant

  R4 = 1/oM;       % R4=2.95e-2;
  R6 = tauR;

  
% swash plate angles in wind axes
%  a1sw= fiC*cos(fiT)*vC/xNN - fiC*sin(fiT)*uC/xNN; % (69) s.17
  a1sw= fiC*(cos(fiT)*vC - sin(fiT)*uC)/xNN;
  
%  b1sw= -fiC*sin(fiT)*vC/xNN - fiC*cos(fiT)*uC/xNN; % (70) s.17
  b1sw= -fiC*(sin(fiT)*vC + cos(fiT)*uC)/xNN;

  % pC: main rotor shaft roll rate in wind control axis 
  pC= -(pB*uC+qB*vC)/xNN;  % (35) s.14
  % qC: main rotor shaft pitch rate in wind control axis 
  qC=(- qB*uC+pB*vC)/xNN;  % (36) s.14

  % ground distance for ground effect calculation
  zB= 5000;
  
  [cTs2, cT2, lamNf, lamI, kG, jItr]= clcViNf(mu, wC, oM, thtaN, pC, zB, cT);
  
  TMRx = cT2*(rhoAir*aDiskMR*(oM*rRot)^2);
  Qmr2=0;
  [aN, a1r, a2r, b1r, b2r, dQ] = feval(clcRtAngl_hndl,mu, lamNf, thtaN, oMR, TMRx);

  % fprintf('Nf cTrl  cTs %10.5f    cT %10.5f   lamNf %10.5f   lamI %10.5f   jItr %3d\n',cTs2,cT2,lamNf, lamI, jItr);
% coning angle aN 
%  aN = eq1(mu, lamNf, thtaN);

% a1 and b1 are rotor flap angles refered to the no feathering axis
  a1= (a1r + (R4*pC - R6*qC))/(1-0.5*mu^2);  % angular terms from eq (6)
  % fprintf('aN %10.2f\n',aN);

  b1= (b1r - (R4*qC + R6*pC))/(1+0.5*mu^2);      % angular terms from eq (7)
  % fprintf('Nf #--# a1 %10.5f    b1 %10.5f\n',a1,b1);
  % fprintf('Nf #--# a1 %10.5f   b1 %10.5f\n',rad2deg(a1),rad2deg(b1));
  

     
  lamNfOut  = lamNf; 
  
  vCt   =[0 0 0];
  omCt  =[0 0 0];
  anglC =[0 0 0 0 0 0 0 0 0 0 0];

  vCt(1)  = uC;  % longitudinal velocity, control axis
  vCt(2)  = vC;  % lateral velocity, control axis
  vCt(3)  = wC;  % vertical velocity, control axis
               
  omCt(1) = pC;  % control axis roll rate
  omCt(2) = qC;  % control axis pitch rate
  omCt(3) = rB;
 
  anglC(1)  = a1;
  anglC(2)  = A1S;
  anglC(3)  = b1;
  anglC(4)  = B1S;
  anglC(5)  = thtaN;
  anglC(6)  = lamI;
  anglC(7)  = mu;
  anglC(8)  = oM;
  anglC(9)  = aN;
  anglC(10) = R4;
  anglC(11) = R6;
  anglC(12) = lamNf;
% angle of no feathering (control) axis
  anglC(13) = -atan(wC/xNN);
  anglC(14) = dQ;
  anglC(15) = cT;
  
% fprintf('%10.5f   %10.5f    %10.5f\n',vCt(1),vCt(2),vCt(3));
% fprintf('%10.5f   %10.5f    %10.5f\n',omCt(1),omCt(2),omCt(3));
% fprintf('%10.5f   %10.5f    %10.5f\n',anglC(1),anglC(2),anglC(3),anglC(4),anglC(5));
 
  
