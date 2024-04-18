function out = cQit(lamNf, mu, thtaN, DvLp, DvLc, PvL, oM)
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
% this modul performs the inflow calculation using Padfields
% iterative scheme
%
% Input:
% lamNf = velocity component perpendicular to no feathering 
%         (control) plane devided by tip speed oMR
% mu    = advance ratio  -> V*cos(alfaD)/oMR
% thtaN = collective pitch at blade root
% DvLp  = drag (except rotor) devided by lift
% DvLc  = aircraft weight times tangent of climb anlge devied by lift
% PvL   = rotor power devided by lift
% oM    = angular rotor speed
%
% Output:
% cQout   = torque/power equilibrium value
%

  global  aLift sigma aDiskMR thtaClmb rhoAir rRot mHeliG...

cT2sa = eq6(mu, lamNf, thtaN);
mucT2sa = mu*cT2sa;

cT = 0.5*sigma*aLift*cT2sa;
DvLi= cT/(2*mu*sqrt(lamNf*lamNf+mu*mu)); % 


tanAlfaNf= lamNf/mu + DvLi;

alNf = atan(tanAlfaNf);

cL = 2*cT*cos(alNf)^3/mu;

DvLi2 = 0.25*cL;

DvLN = eq13(mu,lamNf,thtaN)/mucT2sa;

xh = DvLN + DvLi;

DvLNq = (eq9a(mu,lamNf,thtaN) - eq11d(mu,lamNf,thtaN));

xh = rhoAir*aDiskMR*(oM*rRot)^2*rRot*DvLNq/mHeliG;

% out2 = (DvLp-DvLc)*aLift*mucT2sa + DvLi + DvLN - PvL;

cQout  = (DvLp-DvLc)*aLift*mucT2sa + DvLi + DvLNq - PvL;

