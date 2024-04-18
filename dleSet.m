function dleSet(outP)
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
% This module sets the global variables for longitudinal 
% cyclic to horizontal stabilizer angle gearing curve
% 
% Input:
% outP=1 screen output for test purposes
% 
% Output
% global variables for longitudinal cyclic to horizontal 
% stabilizer angle gearing curve
%
% 
%
global dle2dlsP dleMax dleMin fiHsMin fiHsMax iRigHS
%
% data points for the longitudinal cyclic to horizontal stabilizer angle 
% gearing curve. Minus signs have been added for the cm values
%
outP=0;
	dLngCm(1) =-16.38;   dLngInch(1) =-6.45;    deHS(1) = 0.0224;
	dLngCm(2) =-15.25;   dLngInch(2) =-6.00;    deHS(2) = 0.0174;
	dLngCm(3) =-12.70;   dLngInch(3) =-5.00;    deHS(3) = 0.0000;
	dLngCm(4) =-10.16;   dLngInch(4) =-4.00;    deHS(4) =-0.0192;
	dLngCm(5) =-7.620;   dLngInch(5) =-3.00;    deHS(5) =-0.0384;
	dLngCm(6) =-5.080;   dLngInch(6) =-2.00;    deHS(6) =-0.0541;
	dLngCm(7) =-2.540;   dLngInch(7) =-1.00;    deHS(7) =-0.0690;
	dLngCm(8) = 0.000;   dLngInch(8) = 0.00;    deHS(8) =-0.0820;
	dLngCm(9) = 2.540;   dLngInch(9) = 1.00;    deHS(9) =-0.0850;
	dLngCm(10)= 5.080;   dLngInch(10)= 2.00;    deHS(10)=-0.0803;
	dLngCm(11)= 7.620;   dLngInch(11)= 3.00;    deHS(11)=-0.0628;
	dLngCm(12)= 10.16;   dLngInch(12)= 4.00;    deHS(12)=-0.0300;
	dLngCm(13)= 12.70;   dLngInch(13)= 5.00;    deHS(13)= 0.0035;
	dLngCm(14)= 15.25;   dLngInch(14)= 6.00;    deHS(14)= 0.0593;
	dLngCm(15)= 16.38;   dLngInch(15)= 6.45;    deHS(15)= 0.0942;
	
	dleMin=dLngCm( 1); fiHsMin=deHS(1);
	dleMax=dLngCm(15); fiHsMax=deHS(15);

	nPoly=3;
	nPolyDle=nPoly;
	
    
    % to put the  H-stab gearing curve into effect the H-stab rigging angle
    % has to be -1.5e-6. 
    
  if abs(1.5e-6 + rad2deg(iRigHS)) < 1.0e-9
	 dle2dlsP = polyfit(dLngCm,deHS,nPoly);
  else
     dle2dlsP = [0 0 0 iRigHS];
  end
	if (outP > 0)
	  nVal=250;
	  dLng=(dLngCm(15)-dLngCm(1))/nVal;
	  x=(dLngCm(1):dLng:dLngCm(15));
	  
	  
		for i=1:15
		  yInt(i) = polyval(dle2dlsP,dLngCm(i));  
		  dlta = (yInt(i)-deHS(i))/deHS(i);  
		  fprintf('i %2d  %10.4f   %10.6f    %8.1f\n',i,deHS(i),dlta,dlta*100)  
		end
		
		figure(1)
		plot(dLngCm,yInt)
		hold on
		
		for i=1:15
			plot(dLngCm(i),deHS(i),'+r');
		end
		
		for i=1:nPoly+1
		   fprintf('%d%s%15.9e\n',i,'  ',dle2dlsP(i));
		end   
	end
	
	 


 









