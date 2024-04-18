function [fMR, mMR, oMDotMR] = MainRotorFrcMmtFunc(vCt, omCt, Qmr, anglC)
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
%         --- Version 1.1 ---
%  rigging angle added to force components 102610
%
% This routine calculates main rotor forces and moments
%
% Input:
% vCt   : velocities in control plane (swash plate) coordinates
% omCt  : angular velocities in control plane (swash plate) coordinates
% Qmr   : moment transmitted from rotor to fuselage due to friction
% anglC: control angles
% anglC(1) = a1   : first harmonic longitudinal flapping coefficient
% anglC(2) = A1Swp  : pilot lateral swash plate input
% anglC(3) = b1   : first harmonic  lateral flapping coefficient
% anglC(4) = B1Swp  : pilot longitudinal swash plate input
% anglC(5) = thtaN: collective pitch angle
% lamNf: inflow calculated using Padfields iterative scheme
%
% Output
% fTR : main rotor forces
% mTR : main rotor moments
% oMDotMr : rotor speed derivative
%
%
%
  global rhoAir sigma aLift aDiskMR...
  eRot rRot xGBlade wBlade iBlade nB iRigMR...
  hMR xCG kG dDrgN aN thta1 gloVars cDrgGlo firstOutMR

% R1 R2 R3 R4 R5 R8 R9  
  uC    = vCt(1);% longitudinal velocity, Body axis
  vC    = vCt(2);% lateral velocity, Body axis
  wC    = vCt(3);% vertical velocity,Body axis

  pC    = omCt(1);% body axis roll rate
  qC    = omCt(2);% body axis pitch rate
  rC    = omCt(3);% body axis yaw rate

  a1    = anglC(1);
  A1Swp = anglC(2);
  b1    = anglC(3);
  B1Swp = anglC(4);
  thtaN = anglC(5);

  mu    = anglC(7);
  oM    = anglC(8);
  aN    = anglC(9);
 
  R4    = anglC(10);
  R6    = anglC(11);
  lamNf = anglC(12);
  alfaNf = anglC(13);
  dQ     = anglC(14);
  
  
  oMR = oM*rRot;
  RN  = rhoAir*pi*rRot^2*(oMR)^2;
  R1  = 0.5*sigma*aLift*RN; % R1=17.1e5

  R2=1/(4*oM);   % R2=7.37e-3
  R3=0.5/aLift;  % R3=0.087;


  % thta75 is collective pitch angle at 75% radius as proposed by Bramwell
  thta75 = thtaN + 0.75*thta1;
  xNN=sqrt(uC*uC+ vC*vC);

  
  % in report 73254 R1*xh equals R1*(2*cT/(sigma*aLift)) = cT2sA,
  % thus xh=2*cT/(sigma*aLift)
  xh  = eq6(mu,lamNf,thtaN);
  cTs = 0.5*aLift*xh;
  Tmr = R1*(xh +  R2*mu*pC);
  %
  Hmr =(R3*dDrgN*mu + a1*(thta75/3 +0.75*lamNf*(1-3/8*mu^2) + 0.25*mu*a1) - 0.5*mu*lamNf*thta75*((1-2*mu)/(3*pi)) -0.25*lamNf^2*mu...
       -3/8*lamNf*mu^2*(1+3.333*mu) - aN*(b1/6 - 0.25*aN*mu) - 0.1*R4*qC*(aN/6 + mu*b1/16)...
       -R4*pC*(thta75/6 + 0.5*lamNf + mu*a1/16));
   
  Hmr = 1.15*R1*Hmr; 

  
  DvLN= aLift/(mu*2*cTs)*eq13(mu,lamNf,thtaN); % cos(alfaNf)*
  DvLi= 0.5*cTs*sigma/(mu*sqrt(mu*mu + lamNf*lamNf));
  Hmr2 = Tmr*(DvLN + DvLi);
  
  Ymr = aN*(a1*(1/6-mu^2) - 1.5*mu*(lamNf+0.5*thta75)) + b1*(thta75/3*(1+1.5*mu^2) + 0.75*lamNf + 0.25*mu*a1); % ...
  Ymr = R1*(Ymr + R4*(qC*(thta75/6 +0.5*lamNf + 7/16*mu*a1) - pC*(aN/6 -5/16*mu*b1)));

% With the aircraft moving in the positive x direction the airspeed
% past the hub uc is in the negative x direction for the system of
% coordinates of airspeed used in 73254. Thus uC is < 0 and so is
% the component Hmr*uC
  Xc= (Hmr*uC + Ymr*vC)/xNN; % (14) p.11
  Yc= (Hmr*vC - Ymr*uC)/xNN; % (15) p.11
  Zc= -Tmr;  % (16) p.11

  % Zc points in the negative direction thus -Zc is > 0,
  %   cos(thta), 0,-sin(thta);...
  %       0,     1,      0;...
  %   sin(thta), 0, cos(thta)];
    
  Xr=Xc*cos(B1Swp) - Zc*sin(B1Swp);  % (17) p.11
  Yr=Yc*cos(A1Swp) - Zc*sin(A1Swp);  % (18) p.11
  Zr=Xc*sin(B1Swp) + Zc*cos(B1Swp)*cos(A1Swp) + Yc*sin(A1Swp);  % (19) p.11
  
  fMR = RTrfY(-iRigMR,Xr,Yr,Zr)'; %'
  
  if firstOutMR > 0
    gloVars.Fx=fMR(1); gloVars.Fz=fMR(3);
    gloVars.FxDsk = Hmr; gloVars.FzDsk=Tmr;
    gloVars.Fx= Xr; gloVars.Fz=Zr;
    gloVars.lamNf = lamNf;
    gloVars.mu = mu;
    firstOutMR=0;
  end
  
  % hub moment for offset hinge blades. If the term in brackets 
  % is (B1Swp - a1) a negative sign is required 
  cMRe=  0.5*nB*wBlade*xGBlade*eRot*rRot*oMR^2;
  MRe = -cMRe*(B1Swp - a1);
  LRe =  cMRe*(A1Swp - b1);

  % offset hinge hub moment also acts about x axis
  Lr =  fMR(2)*hMR + LRe;

  % Zr is negative since body z coordinate points downward. For xCG > 0 the 
  % Y-moment contribution from Zr is negative therefore the Zr (=fMR(3)) 
  % term below has positive (+) sign
  Mr = -fMR(1)*hMR + fMR(3)*xCG + MRe;
  Nr = -fMR(2)*xCG;

  mMR = [Lr Mr Nr];
  % QMRtrns = -355;
  QMRtrns = 0.0;
  
  
  % QmRes = 0.5*sigma*rhoAir*aDiskMR*oMR^2*rRot*( eq9a(mu,lamNf,thtaN,oM) - eq11d(mu,lamNf,thtaN,oM) ) + Qmr + QMRtrns;
  QmRes = 0.5*sigma*rhoAir*aDiskMR*oMR^2*rRot*( dQ ) + Qmr + QMRtrns;

  
  % rC turns clockwise when viewed from above since the Z axis points downward
  % while the rotor turns counter clockwise
  oMDotMR = QmRes/iBlade/nB - rC;
  