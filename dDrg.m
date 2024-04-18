function [dDrgN,dDrg1,dDrg2]=dDrg;
% 
% Copyright 2010/2011 Juergen Humt
% 
% This file is part of RotaryWingSim 2.0
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
% This modul calculates coefficients for a three term drag polar
% using the method proposed in naca 716
%
% Input:
% none
% 
% Output:
% default coefficients for three term drag polar
%
global thta1 aLift cLftGlo cDrgGlo

kNDrg  =  0.0003;
k1Drg  = -0.0025;
k2Drg  =  0.0229;


cDrgMinStd = 0.007;
cLMaxStd   = 1.74;
cLOpt      = 0.08;
ReStd      = 8.16e6;
Re         = 2.00e6;

% reynolds number correction
cDrgMin    = cDrgMinStd*(ReStd/Re)^0.11; % fprintf('%s%10.4f\n','cDrgMin    ',cDrgMin);
cLMax      = cLMaxStd - 0.27;   % fprintf('%s%10.4f\n','cLMax      ',cLMax);

cD1=(cLMax-cLOpt);


dDrgN = cDrgGlo*(cDrgMin + kNDrg - k1Drg*cLOpt/cD1 + k2Drg*cLOpt^2/cD1^2);  % fprintf('%s%10.4f\n','dDrgN      ',dDrgN);
cDrgGlo2 =1.0;
dDrg1 = cDrgGlo*( k1Drg/cD1 -2*k2Drg*cLOpt/cD1^2)*aLift;  % fprintf('%s%10.4f\n','dDrg1      ',dDrg1);
dDrg2 = cDrgGlo*(k2Drg/cD1^2)*(aLift)^2; % fprintf('%s%10.4f\n','dDrg2      ',dDrg2);

