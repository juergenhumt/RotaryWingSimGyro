function [ERROR_FLAG] = readAirfoilData(dataFileName,jRe)
% 
%
% Copyright 2010/2011 Juergen Humt
% 
% This file is part of rotaryWingSim.
% 
%     rotaryWingSim, is free  software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by the 
%     Free Software Foundation, either version 3 of the License or any later 
%     version.
% 
%     rotaryWingSim is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License along 
%     with rotaryWingSim.  If not, see <http://www.gnu.org/licenses/>.
%
%
%  Version 1.0 of this program was based on report naca-tm-73254 
%
%         --- Version 1.3 ---
%
%
% This file reads airfoil data
%
% Input:  
% dataFileName   (Name of file containing airfoil data)    (string)
%
% Output:
% ERROR_FLAG  (1 if error has occurred 0 otherwise)
% ERROR_FLAG is set to 1 if there is a problem opening airfoil data file
%
% Output to workspace:
% numData     number of data points read for this airfoil
% AlphaZL (Zero lift AoA)       [rad]

% airfoile data file format ("name.pol")
% numRe
% Re
% numData
% AlphaEnd  [degrees]
% AlphaZL   [degrees]
% AlphaData [degrees] ClData CdData CmData 
%
global AlphaEnd AlphaZL numData ClData CdData CmData ReData numRe
ERROR_FLAG = 0;

dataFile = fopen(dataFileName);          % Open data file for chosen airfoil

if ( dataFile ~= -1 ),
    
ReData(1)=0;
AlphaEnd(1)=0;

% skip nine lines in file
readDmy(dataFile,8);
d1 = fscanf(dataFile, '%s', 5);

% read line containing Reynolds number
d1 = fscanf(dataFile, '%s', 3);
ReData(jRe) = sscanf(d1,'%f');

% fprintf('%s%f\n','Re ',ReData(jRe));
readDmy(dataFile,4);

iDat=0;
while feof(dataFile) == 0

  TempData = fscanf(dataFile, '%f', 7);  % Xfoil spits out 7 types of data in its polars
  if 7 == length(TempData)      
    iDat=iDat+1;
    fi = deg2rad(TempData(1));
    ClData(iDat,1) = fi;
    CdData(iDat,1) = fi;
    CmData(iDat,1) = fi;
    
    ClData(iDat,jRe)  = TempData(2);
    if ( abs(ClData(iDat,jRe)) < 1.0e-6 )
        AlphaZL(jRe)= ClData(iDat,jRe);
    end
    CdData(iDat,jRe)  = TempData(3);
    CmData(iDat,jRe)  = TempData(5);
  end 
end


fclose(dataFile);


xtndDat(jRe);
[nDat jh]=size(ClData);
numData(jRe)=nDat ;


reX=1.0e6;
xh=1.05;

for k=3:4
  ReData(k)=reX;
  AlphaEnd(k)=AlphaEnd(jRe);  
  for i=1:nDat
    ClData(i,k) = xh*ClData(i,2);
    CdData(i,k) = xh*CdData(i,2);
    CmData(i,k) = xh*CmData(i,2);
  end
  reX=2.0e6;
  xh=xh*xh;
end

numRe=4;

else
  disp(sprintf('ERROR AIRFOIL.M: Error opening data file %s!', Airfoil));
  Airfoil = '';
  ERROR_FLAG = 1;
end
