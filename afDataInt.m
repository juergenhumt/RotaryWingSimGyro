function [cL, cD, cM, ERROR_FLAG] = afDataInt(alfaIn, ReIn)
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
%         --- Version 2.0 ---
%
% 
% 
%  This modul calculates wing coefficients from polynomial values
%
%  Input:
%  alfaIn: angle of attack of airfoil
%  ReIn  : Reynolds number of flow
%
%  Output:
%  cL: section lift coefficient
%  cD: section drag coefficient
%  cM: section momement coefficient
%  ERROR_FLAG: is 1 if section is stalled
%
global numRe numData AlphaEnd ReData ClData CdData CmData cDrgGlo cLftGlo CdMax


   ERROR_FLAG=0;
   if (ReIn < ReData(numRe))
       if (ReIn < ReData(2))
           reFact=0;
           % reInd must be 3 since interpolation is 
           % perfromed between reInd and reInd-1
           % and reInd=2 is the lowest Re in the table
           % modHu 170911 reInd was 2 before
           reInd=3;
       else 
           j=3; 
           while (ReData(j) < ReIn)
             j=j+1;
           end
           dRe = ReData(j)-ReData(j-1);
           reInd=j;
           reFact= (ReIn-ReData(j-1))/dRe;
       end
   else
       reFact= 1;
       reInd = numRe;
   end    


   if( alfaIn > AlphaEnd(reInd) ),		% Check to see if element is stalled
			% sprintf('AIRFOIL.M: Element %i stalled at AoA = %f!', i, alfaIn)
			% disp('AIRFOIL.M: Using cl=0, cd=1.0, and cm=0.')
			ERROR_FLAG = 1;
			cL = 0;
			cD = CdMax;
			cM = 0;

	else
        numDat = numData(2);
        % fprintf('numData %8.0f\n',numDat);
        
        if( alfaIn < ClData(1,1)),
			% sprintf('AIRFOIL.M: Element %i at AoA = %f is below experimental data range!', i, alfaIn)
			% sprintf('AIRFOIL.M: Using cl, cd, and cm at the minimum AoA, %f degrees', rad2deg(ClData(1,1)))
            cL = ClData(1,reInd-1) + reFact*(ClData(1,reInd)-ClData(1,reInd-1));
            cD = CdData(1,reInd-1) + reFact*(CdData(1,reInd)-CdData(1,reInd-1));
            cM = CmData(1,reInd-1) + reFact*(CmData(1,reInd)-CmData(1,reInd-1));
			ERROR_FLAG = 2;
		% If local AoA is above maxAlpha (but below AlphaEnd), interpolation on Re is impossible, 
        % so just use values from the data set of maxAlpha
		elseif( alfaIn > ClData(numDat,1)),
			% sprintf('AIRFOIL.M: Element %i at AoA = %f is above experimental data range!', i, alfaIn)
			% sprintf('AIRFOIL.M: Using cl and cd at the maximum AoA, %f degrees', rad2Deg(ClData(numDat)))
            cL = ClData(numDat,reInd-1) + reFact*(ClData(numDat,reInd)-ClData(numDat,reInd-1));
            cD = CdData(numDat,reInd-1) + reFact*(CdData(numDat,reInd)-CdData(numDat,reInd-1));
            cM = CmData(numDat,reInd-1) + reFact*(CmData(numDat,reInd)-CmData(numDat,reInd-1));
            ERROR_FLAG = 3;
        else
		%  interpolate now
            alfaInt =alfaIn;
            m    =fIDX(ClData(:,1),alfaInt);
            % alfa =rad2deg(alfaInt);
            
			cL = intDat(ClData,alfaInt,m,reInd,reFact);	% cl
			cD = intDat(CdData,alfaInt,m,reInd,reFact);	% cl
			cM = intDat(CmData,alfaInt,m,reInd,reFact);	% cl
            % fprintf('%s%12.8f\n','cL    ',cL);
            % fprintf('%s%12.8f\n','cD    ',cD);
            % fprintf('%s%12.8f\n','cM    ',cM);
            
        end 

% else,
% disp(sprintf('ERROR AIRFOIL.M: Unknown problem with element %i at AoA = %f!', i, alfaIn));
% ERROR_FLAG = 1;
        cL=cL;
        cD=cDrgGlo*cD;
    end
