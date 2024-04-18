function [A_cf,B_cf,C_cf,D_cf,E_cf]= coeffnc (Xu,Xw,Xq,Xwd,Zu,Zw,Zq,Zwd,Mu,Mw,Mwd,Mq,Mthta)
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
% stability quartic 
% 
% Input:
%  stability derivatives
% 
% Output
%  coefficients A, B, C, D, E of stability quartic
%

global Iyy m gEarth

% s^4
A_cf=1;

% s^3
B_cf=-(Mq +Zw +Xu);

% s^2
C_cf= (-Mu*Xq+Xu*Mq-Zu*Xw+Xu*Zw-Mw*Zq+Zw*Mq);

% s^1
D_cf=(Mu*gEarth+Mu*Zw*Xq-Mu*Xw*Zq+Zu*Xw*Mq-Xu*Zw*Mq+Xu*Mw*Zq-Zu*Mw*Xq);

% s^0
E_cf= (Zu*Mw-Mu*Zw)*gEarth;

% N = -Xu -Zw
% P =  Xu*Zw - Xw*Zu
% Q = -(mu*cos(alfaD) + Zq)*Xu - cTs*sin(gam) + Zu*Xq
% R = -cTs*( Zw*cos(gam) - Xw*sin(gam))
% S =  cTs*cos(gam) - Xw*(mu/cos(alfaD) + Zq) - Zw*Xq 
% T = -cTs*(Zw*cos(gam) - Xw*sin(gam))


% w = -Mw
% nu= -Mq
% Xi= -Mwdot
% H = -Mu

% A = 1
% B = N + nu + Xi*(mu/cos(alfaD) + Zq)
% C = P + nu*N + Xi*() + w*(mu/cos(alfaD) - Zq) + H*Xq
% D = nu*P + Xi*R * Q*w - S*H
% E = w*R - H*T


% B = (-Xu -Zw) + -Mq + -Mwdot*(mu/cos(alfaD) + Zq)
% C = Xu*Zw - Xw*Zu + -Mq*N + Xi*(-(mu*cos(alfaD) + Zq)*Xu - cTs*sin(gam) + Zu*Xq) + w*(mu/cos(alfaD) - Zq) + H*Xq
% D = -Mq*(Xu*Zw - Xw*Zu) + Xi*R + (-(mu*cos(alfaD) + Zq)*Xu - cTs*sin(gam) + Zu*Xq)*w - S*H
% E = w*R - H*T


