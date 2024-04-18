function zetaSet(outP)
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
% This module calculates the polynomial coefficents for the two curves for 
% induced velocity  amplification factor given by Bramwell
%
% Input
%   outP   if outP greater zero then generate plot of zeta curves
% 
% Output 
%   global variables for polynomial coefficients
% 

global zPA zPB 
nPoly =6;
xCDf=...
 [-0.31
 -0.25
 -0.2
 -0.15
 -0.1
 -0.05
 -0.0
  0.05
  0.1
  0.15
  0.2
  0.25
  0.3];


yCDf =...
  [1.05
   1.34
   1.66
   1.98
   2.1
   2.05
   1.95
   1.82
   1.67
   1.53
   1.39
   1.24  
   1.14];

% fprintf('%s\n\n','### polyfit ###');
zPA = polyfit(xCDf,yCDf,nPoly);

if outP > 0
close all
for i=1:13
  yInt(i) = polyval(zPA,xCDf(i));  
  dlta = (yInt(i)-yCDf(i))/yCDf(i);  
  fprintf('%s%10.0f%s%10.2f%s%10.5f\n','d  ',i,'  ',xCDf(i),'  ',dlta)  
end

figure(1)
plot(xCDf,yInt),axis([-0.3,0.3,0,3.0])
hold on

for i=1:13
	plot(xCDf(i),yCDf(i),'+r');
end

for i=1:nPoly+1
   fprintf('%d%s%15.9g\n',i,'  ',zPA(i));
end   
end

xCDf= [-0.31
-0.25
-0.2
-0.15
-0.1
-0.05
-0.0
 0.05
 0.1
 0.15
 0.2
 0.25
 0.3];


yCDf = [1.3
   1.41
   1.62
   1.70
   1.70
   1.65
   1.50
   1.42
   1.34
   1.26
   1.18
   1.10  
   1.06];

% fprintf('%s\n\n','### polyfit ###');
zPB = polyfit(xCDf,yCDf,nPoly);

if outP > 0
for i=1:13
  yInt(i) = polyval(zPB,xCDf(i));  
  dlta = (yInt(i)-yCDf(i))/yCDf(i);  
  fprintf('%s%10.0f%s%10.2f%s%10.5f\n','d  ',i,'  ',xCDf(i),'  ',dlta)  
end

plot(xCDf,yInt,'-g')
for i=1:13
	plot(xCDf(i),yCDf(i),'xr');
end

for i=1:nPoly+1
   fprintf('%d%s%15.9g\n',i,'  ',zPB(i));
end   
end
