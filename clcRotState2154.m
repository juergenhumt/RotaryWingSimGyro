function [thtaN, thtaNDeg, lamNf, gamDeg, alfaNf, alfaNfDeg, alfaDeg, vFlight, vH, vV, nRotMin, oMR, Tmr]  = clcRotState2154(cTs, mu, cDj)
%
% RotaryWingSimGyro
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
% This routine calculates the rotor inflow for autorotation and the collective
% pitch requried for the given cTs. From there several other paramters are 
% calculated. Rotor parameter include a tip jet engine, the formulae are from
% naca report 2154 "Autorotative Performance of a Tip Jet Helicopter"
%
%
global bTl gamRot aLift thta1 cDrgGlo cLftGlo cHMR sigma rhoAir...
  dDrgN dDrg1 dDrg2 thta1 rRot aDiskMR fpA mHeliG...
  dDrgJet dJet aJet rOutr rInnr rCntr...
  t31 t32 t33...
  t41 t42 t43 t44 t45 t46 t47 t48 t49 t41N...
  t51 t52 t53 t54 t55 t56 t57 t58 t59 t51N

  fpA2154 = fpA/aDiskMR;

     if cDj > 0
      xh=(0.5*mu^2*(rOutr - rInnr) + (rOutr - (rInnr^3/rOutr^2))/3);
      rJR = (rCntr/rOutr);
      xh=xh/(rJR^2+0.5*mu^2);

      dcDj = cDj - dDrgJet*cHMR/aJet*xh ;

      cQj = dcDj *aJet/aDiskMR*rJR*(rJR^2 + 0.5*mu^2);
      DvLj= aJet/aDiskMR*dcDj/(mu*cT)*(rJR^3+1.5*mu^2*rJR);
    else
      cQj = 0;
      DvLj= 0;
    end

 
  
    coeff716(mu);

    aLiftLoc=aLift;

    t3a= t32/t31;
    cTp= (2*cTs/aLiftLoc -t33*thta1)/t31;
    cT = cTs*sigma;

  % NACA 716 s 217
    f2 =  t41*t3a^2 - t42*t3a + t44;
    f1 = -2*t41*t3a*cTp + t42*cTp - t43*t3a*thta1 + t45*thta1;
    fN =  t41*cTp^2 + t43*cTp*thta1 + t46*thta1^2;

    f2 = aLiftLoc*f2;
    f1 = aLiftLoc*f1;
    fN = aLiftLoc*fN;


    g2 =  dDrg2*(-t56*t3a + t58 + t55*t3a^2);
    g1 =  dDrg1*(-t3a*t52 + t53) + dDrg2*(t56*cTp - t57*t3a*thta1 + t59*thta1 -2*cTp*t55*t3a);
    gN =  dDrgN*t51 + dDrg1*(t52*cTp+t54*thta1)...
         + dDrg2*(t55*cTp^2+t57*cTp*thta1+t51N*thta1^2) + 2*cQj/sigma;


    pE = (g1-f1)/(g2-f2);
    qE = (gN-fN)/(g2-f2);


  % lamda page 217 naca716
    thtaN = -pE/2+sqrt(pE^2/4-qE);

    thtaNDeg= rad2deg(thtaN);
   % fprintf('thta %9.6f°\n',thtaNdeg);

    lamNf= (2*cTs/aLift - t32*thtaN - t33*thta1)/t31;


    % disk angle
    DvLi=cT/(2*mu*sqrt(lamNf*lamNf+mu*mu));
    tanAlNf=lamNf/mu + DvLi;


    alfaNf = atan(tanAlNf);
    alfaNfDeg=rad2deg(alfaNf);

    a1 = eq2(mu,lamNf,thtaN);

    alfaD = alfaNf + a1;
    alfaDeg = rad2deg(alfaD);

    DvLj=DvLj*cos(alfaNf);
    DvLN= cos(alfaNf)*aLift/(mu*2*cTs)*eq13(lamNf,thtaN);
    
    DvLi = DvLi*cos(alfaNf);
    DvLp = fpA2154*mu^2/(2*cTs*sigma*cos(alfaNf)^2);

    DvLg = DvLN + DvLi + DvLp + DvLj;


    siGam= DvLg/sqrt(1 - DvLp^2 + 2*DvLp*DvLg);


    gamma = asin(siGam);
    gamDeg = rad2deg(gamma);

    Tmr = mHeliG*tan(gamma)/DvLg;

    oMR = sqrt(Tmr/(cTs*sigma*aDiskMR*rhoAir));

    nRotMin = oMR/(2*pi*rRot)*60;

    vFlight= mu*oMR/cos(alfaNf);

    vH= vFlight*cos(gamma);
    vV= vFlight*siGam;


    k=-1;
    if k > 0
    cT2sA = t31*lamNf+t32*thtaN+t33*thta1;

    cTs = 0.5*cT2sA*aLiftLoc;
    cT  = cTs*sigma;

    end
