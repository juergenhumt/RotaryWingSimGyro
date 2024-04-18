function xtndDat(kRe)
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
% This module calculates extensions for the wing profile Cl, Cd and Cm curves
%
% Input
%   kRe  number of data set, each data file contains data for one Reynolds number
% 
global AlphaEnd AlphaZL numData ClData CdData CmData CdMax ReData numRe

[numCl kCl]=size(ClData);

kRe=2;
% fprintf('cl  %10.5f  %10.5f\n', ClData(1,1), ClData(1,kRe));
% fprintf('cl  %10.5f  %10.5f\n', ClData(numCl,1), ClData(numCl,kRe));

dyl= ClData(numCl,kRe) - ClData(numCl-1,kRe);
dyd= 0.15*CdData(numCl,kRe);


jCl = numCl;


xhl =  ClData(numCl,1) - AlphaZL(kRe);
dxl = xhl/10;



if (dyl > 0)
% if cl is not decreasing at the end dummy points are added
  dyl = dyl/50;
  jCl=jCl+1;
  ClData(jCl,1) = ClData(numCl,1)+dxl;
  ClData(jCl,kRe) = ClData(numCl,kRe)+dyl;
  
  
  CdData(jCl,1) = ClData(numCl,1)+dxl;
  CdData(jCl,kRe) = CdData(numCl,kRe)+dyd;
  
  CmData(jCl,1) = ClData(numCl,1)+dxl;
  CmData(jCl,kRe) = 0.8*CmData(numCl,kRe);
  
  
  jCl=jCl+1;
  ClData(jCl,1) = ClData(numCl,1)+ 2*dxl;
  ClData(jCl,kRe) = ClData(numCl,kRe);
  
  CdData(jCl,1) = ClData(numCl,1)+2*dxl;
  CdData(jCl,kRe) = CdData(numCl,kRe)+4*dyd;
  
  CmData(jCl,1) = ClData(numCl,1)+2*dxl;
  CmData(jCl,kRe) = 0.5*CmData(numCl,kRe);

  jX=3;
else
  jX=1;
end

jCl=jCl+1;
ClData(jCl,1) = ClData(numCl,1)+ jX*dxl;
ClData(jCl,kRe) = 0.01*ClData(numCl,kRe);

CdData(jCl,1) = ClData(numCl,1)+jX*dxl;
CdData(jCl,kRe) = 0.7*CdMax;

CmData(jCl,1) = ClData(numCl,1)+jX*dxl;
CmData(jCl,kRe) = 0.01*CmData(numCl,kRe);

jC1=jCl;
jCl=jCl+1;
xh=1.3;

ClData(jCl,1) = xh*ClData(jC1,1);
ClData(jCl,kRe) = 0.0;

CdData(jCl,1) = xh*ClData(jC1,1);
CdData(jCl,kRe) = CdMax;

CmData(jCl,1) = xh*ClData(jC1,1);
CmData(jCl,kRe) = 0.0;


AlphaEnd(kRe)=ClData(jCl,1);
[numCl kCl]=size(CdData);
[numCl kCl]=size(CmData);

