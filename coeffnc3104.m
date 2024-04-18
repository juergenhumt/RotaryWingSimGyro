function [A_cf,B_cf,C_cf,D_cf,E_cf]= coeffnc3104 (Xu,Xw,Xq,Xwd,Zu,Zw,Zq,Zwd,Mu,Mw,Mq,Mwd,Mthta,mu,alfaD,gam,cTs)
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
% This module calculates the coefficients for the
% stability quartic based on RaM 3104
% 
% Input:
%  stability derivatives
% 
% Output
%  coefficients A, B, C, D, E of stability quartic
%

global mHeliG rhoAir gEarth sigma aDiskMR rRot

% s^4
% A_cf=1;

% s^3
% B_cf=-(Mq +Zw +Xu);

% s^2
% C_cf= (-Mu*Xq+Xu*Mq-Zu*Xw+Xu*Zw-Mw*Zq+Zw*Mq);

% s^1
% D_cf=(Mu*gEarth+Mu*Zw*Xq-Mu*Xw*Zq+Zu*Xw*Mq-Xu*Zw*Mq+Xu*Mw*Zq-Zu*Mw*Xq);

% s^0
% E_cf= (Zu*Mw-Mu*Zw)*gEarth;


% The forumula for the quatic coefficients in RaM 3104 uses
% several auxiliary terms:
% 
% N = -Xu -Zw
% P =  Xu*Zw - Xw*Zu
% Q = -(mu*cos(alfaD) + Zq)*Xu - cTs*sin(gam) + Zu*Xq
% R = -cTs*( Zu*cos(gam) - Xu*sin(gam))
% S =  cTs*cos(gam) - Xw*(mu/cos(alfaD) + Zq) - Zw*Xq 
% T = -cTs*(Zw*cos(gam) - Xw*sin(gam))


% w = -Mw
% nu= -Mq
% Xi= -Mwd
% H = -Mu


% B = N + nu + Xi*(mu/cos(alfaD) + Zq)
% C = P + nu*N + Xi*() + w*(mu/cos(alfaD) - Zq) + H*Xq
% D = nu*P + Xi*R * Q*w - S*H
% E = w*R - H*T


% B = (-Xu -Zw) + (-Mq + -Mwd*(mu/cos(alfaD)) + Zq)
% C =  Xu*Zw - Xw*Zu + (-Mq*(-Xu -Zw)) + (-Mwd*(-(mu*cos(alfaD) + Zq))*Xu - cTs*sin(gam) + Zu*Xq) + (-Mw*(mu/cos(alfaD)) - Zq) + (-Mu*Xq)
% D = -Mq*(Xu*Zw - Xw*Zu) + (-Mwd*(-cTs*( Zw*cos(gam) - Xw*sin(gam)))) + (-(mu*cos(alfaD) + Zq)*Xu - cTs*sin(gam) + Zu*Xq)*(-Mw)  - (-Mu*(cTs*cos(gam)) - Xw*(mu/cos(alfaD) + Zq) - Zw*Xq)
% E = -Mw*(-cTs*( Zw*cos(gam) - Xw*sin(gam))) - -Mu*(-cTs*(Zw*cos(gam) - Xw*sin(gam)))


mu2 = mHeliG/(rhoAir*gEarth*sigma*aDiskMR*rRot);
cMu  = 0; % Uo has already been added to Zq outside this routine

A_cf =      1;
B_cf = -Xu -Zw -Mq -Mwd*(mu/cos(alfaD) + Zq);
%      +Xu*Zw - Zu*Xw  +Xu*Mq +Zw*Mq       -Mw*Zq-Mw*Un       -Mu*Xq+Mu*w)
C_cf =  Xu*Zw - Xw*Zu -Mq*(-Xu -Zw) -Mw*(cMu*mu/cos(alfaD) + Zq) -Mu*Xq       -Mwd*(- cTs*sin(gam) -(cMu*mu*cos(alfaD) + Zq)*Xu);
%                                                 + Zu*Xq  
%       +Zu*Xw*Mq -Xu*Zw*Mq  +Mu*gE*csG              +Mw*gE*snG     +Mu*Xq*Zw   -Mu*Xw*Un      -Mu*Xw*Zq            +Xu*Mw*Un+        +Xu*Mw*Zq    +Zu*Mw*w-Zu*Mw*Xq 
D_cf =  +Mq*Xw*Zu -Mq*Xu*Zw +Mu*gEarth*cos(gam) +Mw*gEarth*sin(gam) +Mu*Zw*Xq -Mu*(+Xw*(cMu*mu/cos(alfaD) + Zq)) +Mw*Xu*(cMu*mu*cos(alfaD) + Zq)        -Mw*Zu*Xq     -Mwd*(-gEarth*( Zw*cos(gam) - Xw*sin(gam)));

% cTs has been replaced by gEarth since the derivatives are not non-demensional
% E_cf = -Mw*(-cTs*( Zu*cos(gam) - Xu*sin(gam))) + Mu*(-cTs*(Zw*cos(gam) - Xw*sin(gam)));
%        +Zu*Mw*gE*csG-Xu*Mw*gE*snG               + Mu*(Xw*gE*snG-Zw*gE*csG) 
E_cf = -Mw*(-gEarth*( Zu*cos(gam) - Xu*sin(gam))) + Mu*(-gEarth*(Zw*cos(gam) + Xw*sin(gam)));


% test output of symbolic program calculating Stengel's stability matrix
% 
% w=0;  Un=0;
% B_cf= (-Mq-Zw-Xu);
% C_cf= +(Zw*Mq-Mw*Zq-Mw*Un*cMu+Xu*Mq+Xu*Zw-Zu*Xw-Mu*Xq+Mu*w);
% D_cf= +(-Xu*Zw*Mq+Zu*Xw*Mq-Zu*Mw*Xq+Xu*Mw*Un*cMu+Mw*gEarth*sin(gam)-Mu*Xw*Zq-Mu*Xw*Un*cMu+Zu*Mw*w+Mu*gEarth*cos(gam)+Mu*Xq*Zw-Mu*w*Zw+Xu*Mw*Zq);
% E_cf= -Xu*Mw*gEarth*sin(gam)-Mu*Zw*gEarth*cos(gam)+Mu*Xw*gEarth*sin(gam)+Zu*Mw*gEarth*cos(gam);
%  
% A=1;



% +(Xu*Mw*Un
% -Mu*Xw*Un
% +Zu*Mw*w
% -Mu*Zw*w


% coefficients checked
%             772   Brw  0047
% +Zu*Xw*Mq    v     v    v   
% -Xu*Zw*Mq    v     v    v   
% +Xu*Mw*Zq    v     v    v
% -Zu*Mw*Xq    v     v    v
% +Mu*Zw*Xq    -     v    v
% -Mu*Xw*Zq    -     v    v

% +Zu*Mw*gEarth-Mu*Zw*gEarth

