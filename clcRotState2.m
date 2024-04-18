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
% This routine calculatest rotor inflow for autorotation and several
% output parameters starting from that value
%
% Input:
%   thtaN: collective blade angle
%   mu   : adavance ratio
%
% Output:
%   cTs  : thrust coefficient devided by solidity
%   cT   : thrust coefficient
%   lamNf: inflow wrt swashplate (no feathering axis)
%   vi   : induced velocity
%   alfaD: disk angle of attack (rad)
%   alfaDeg: disk angle of attack (deg)
%
%
function [cTs, cT, lamNf, vi, alfaD, alfaDeg, alfaNf]  = clcRotState2(thtaN, mu)
global bTl gamRot thta1 cLftGlo sigma...
  dDrgN dDrg1 dDrg2 thta1 aLift...
  t31 t32 t33...
  t41 t42 t43 t44 t45 t46 t47 t48 t49 t41N...
  t51 t52 t53 t54 t55 t56 t57 t58 t59 t51N

  coeff716(mu);

  % NACA 716 s 217
    f2 =t41;
    f1 =(t42*thtaN+t43*thta1);
    fN =(t44*thtaN^2+t45*thtaN*thta1+t46*thta1^2);

    f2 = aLift*f2;
    f1 = aLift*f1;
    fN = aLift*fN;


    g2 = dDrg2*t55;
    g1 = (dDrg1*t52 + dDrg2*(t56*thtaN+t57*thta1));
    gN =  dDrgN*t51 + dDrg1*(t53*thtaN+t54*thta1)+dDrg2*(t58*thtaN^2+t59*thtaN*thta1+t51N*thta1^2);

    pE = (g1-f1)/(g2-f2);
    qE = (gN-fN)/(g2-f2);


  % lamda page 217 naca716
    lamNf = -pE/2+sqrt(pE^2/4-qE);


    if abs(imag(lamNf)) > 1.0e-10
      fprintf('lamNf %10.5e\',lamNf);
    end


    cT2sA = eq6(mu, lamNf, thtaN);

    cTs   = 0.5*cT2sA*aLift;
    cT  = cTs*sigma;

%     thta75= thtaN+0.75*thta1;
%     cTsBr = 0.25*aLift*(2/3*thta75*(1 + 1.5*mu*mu) + lamNf);

  
  % disk angle
    DvLi=cT/(2*mu*sqrt(lamNf*lamNf+mu*mu)); % is equal to lamI, the induced inflow parameter
    tanAlNf =lamNf/mu + DvLi;
    vi = DvLi*mu;

    alfaNf = atan(tanAlNf);
    
    a1 = eq2(mu,lamNf,thtaN);
    alfaD = alfaNf + a1;
    alfaDeg = rad2deg(alfaD);
