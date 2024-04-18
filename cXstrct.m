function cX = cXstrct(inpFile)
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
% This modul reads coefficients of a polynomial, left and right 
% interval boundaries and an amplification factor. Prescribed 
% values for linearly extrapolating the curves are read as well,
% if present. For other curves for which no prescribed values for 
% linearlinear extraplation are present the extrapolation terms 
% are calculated in this routine. Since the file pointer is not
% altered consecutive reads will read the next block of data 
% from the input file. The block length is the first value after
% the comment line
%
% Input:
% inpFile: file pointer 
% 
% Output:
% cX: data structure
%
zCff(1)=0;


% dummy reading the comment line
  sAux = fgetl(inpFile);
% read prescribed values of linear term outside 
% of interval where coefficients are valid
  sAux = fgetl(inpFile);
  a1  = sscanf(sAux,'%e');
  sAux = fgetl(inpFile);
  a2  = sscanf(sAux,'%e');

% read number of coefficients to follow
  sAux = fgetl(inpFile);
  nDat  = sscanf(sAux,'%e');

  for i=1:nDat+1
    sAux = fgetl(inpFile);
    ret  = sscanf(sAux,'%e');
    % fprintf('%12.8e\n',ret);
    zCff(i) = ret;
  end % for
  
  sAux = fgetl(inpFile);
  x1a  = sscanf(sAux,'%e');
  sAux = fgetl(inpFile);
  x2a  = sscanf(sAux,'%e');

  sAux = fgetl(inpFile);
  cAmp  = sscanf(sAux,'%e');
  
  sAux = fgetl(inpFile);
  cDat  = sscanf(sAux,'%e');
  
  cX.cff = zCff;
  cX.x1  = x1a;
  cX.x2  = x2a;
  cX.nDat= nDat;
  cX.cAmp= cAmp;
  

  x1=cX.x1;
  y1= evlPly(x1, cX);
  
  if (a1==0)
    x2=0.99*cX.x1;
    y2= evlPly(x2, cX);

    a = (y2-y1)/(x2-x1);
    b = (x2*y1 - x1*y2)/(x2-x1);
    
    cX.a1 = a;    
    cX.b1 = b;
  
  else
    cX.a1 = a1;
    cX.b1 = y1 - a1*x1;
  end
  
  
  % ok
  x2=cX.x2;
  y2= evlPly(x2, cX);

  if (a2==0)
    x1=0.99*cX.x2;
    y1= evlPly(x1, cX);
    
    a= (y2-y1)/(x2-x1);
    b= (x2*y1 - x1*y2)/(x2-x1);
    
    cX.a2 = a;    
    cX.b2 = b;
  else
    cX.a2 = a2;
    cX.b2 = y2 - a2*x2;
  end