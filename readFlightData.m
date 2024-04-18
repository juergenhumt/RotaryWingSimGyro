function [thtaN, nRot, oM, oMR] = readFlightData(dataDir,fileName,str2find,outFile)
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
% This module reads data from a named section in an input file the 
% section name is stored in input variable str2Find. Information 
% regarding success of the read operation is written to output file. 
% Some selected values are returned but are returned but these are
% currently not used
%
% Input
%  datadir   name/path of data directory
%  fileName  input file name
%  jRead     number of lines to be skipped
% 
% Ouput
%  thtaN     collective pitch angle
%  nRot      rotor speed
%  oM        angular rotor velocity
%  oMR       tip speed, i.e. angular rotor velocity times rotor radius
% 
global pE vE qE aE usephi nRot nRotMin oM oMR iRigHS...
       rRot xcg lmd mHeli mHeliG thtaClmb PA addLd...
       engineRPM engineHP thtaN A1S B1S qMR...
       f2m i2m lbForce lbMass rhoAir gEarth gloVars rfrncVal itrType...

inFileName= [dataDir fileName '.flt'];
fprintf(outFile,'reading %s from %s\n',str2find,inFileName);

% fprintf('##     Reading ','\n\n%s\n',str2find);
eval(['disp(''# Data   ',inFileName,'   #'')'])

jCount=0;
inpFile = fopen(inFileName,'r');
readCont=1;
if inpFile > 0
  while ((feof(inpFile) == 0) & readCont)
    sAux = fgetl(inpFile);
    jCount=jCount+1;
    % fprintf('%d    %s  \n',jCount,sAux);
    if length(sAux) > 0
    token=strtok(sAux);
    if '%'==token(1)
      k=findstr(str2find,sAux);
      if k > 0
      tmpStr = ['Start reading ' str2find ' section\n'];
      fprintf(tmpStr);
      
      if  strncmpi('Flight Data',str2find,6) > 0

          readDmy(inpFile,1)       %
          readDmy(inpFile,1) % Start position, global coordiantes
          pE(1) =  readLn(inpFile); % pE(1)
          pE(2) =  readLn(inpFile); % pE(2)
          pE(3) =  readLn(inpFile); % pE(3)
          readDmy(inpFile,1) % Initial velocites
          [vAux sAux]  =readLn2(inpFile);  % vE(1)   [knot] -?-
          if (strcmp('knot',sAux) | strcmp('kts',sAux))
            spdCnvFact= 1.805/3.6;
          else
            spdCnvFact= 1.609/3.6;              
          end
          vE(1)= vAux*spdCnvFact;
          
          vE(2) = readLn(inpFile)*spdCnvFact; % 
          vE(3) = readLn(inpFile)*f2m/60; % 
          readDmy(inpFile,1) % Initital rotational speeds
          qE(1) =readLn(inpFile); % qE(1)
          qE(2) =readLn(inpFile); % qE(2)
          qE(3) =readLn(inpFile); % qE(3)
          readDmy(inpFile,1) % Initital orientation
          aE(1) = deg2rad(readLn(inpFile)); % phi
          aE(2) = deg2rad(readLn(inpFile)); % theta
          aE(3) = deg2rad(readLn(inpFile)); % psi

          readDmy(inpFile,1)         
          % additional data
          addLd   =readLn(inpFile)*lbMass;
          mHeli   =mHeli + addLd;
          mHeliG  =mHeli*gEarth;
          
          Temp=readLn(inpFile);           % temperature        [degs F]
          PA=readLn(inpFile);             % pressure  altitude [ft    ]
          if (PA > 0) 
            PA = PA*f2m; 
          end
          if (PA < 0)
            rhoAir = -PA*rhoAir;
          end;
          qMR = readLn(inpFile)*lbForce*f2m;   % main rotor torque [ft-lb   ]

          [vAux sAux]  =readLn2(inpFile);  % rotational  velocity    [rads/sec] -?-
          if vAux > 0
            nRot=0; vTip=0; oM=0;
            if (strcmp('rads/sec',sAux) | strcmp('oM',sAux))
              oM   =vAux;
              nRot = oM/2/pi;
              vTip = rRot*oM;
            end
            if (strcmp('revs/min',sAux) | strcmp('rpm',sAux))
              nRot = vAux/60;
              oM = 2*pi*nRot;
              vTip = rRot*oM;
            end
            if strcmp('fps',sAux)
              vTip = vAux*f2m;
              oM   = vTip/rRot;
              nRot = oM/2/pi;
            end
            
            if ( oM+nRot+vTip < 1.0e-8 )
              fprintf('Input error: omega or nRot or vTip needed\n');
              error(' *** END OF PROGRAM  ***')
            end
  
            oMR = oM*rRot;
            nRotMin = nRot*60;
          end 
          xh    = readLn(inpFile);
          % thtaN = deg2rad(xh);
          xh  = readLn(inpFile);
          xcg = xcg - xh*i2m;  % CG offset                    [ft    ]
          
          % read stability derivative reference values
          readDmy(inpFile,3) 
          k=1;
          rfrncVal=zeros([1 10]);
          while ((feof(inpFile) == 0) & k < 11)
            rfrncVal(k) = readLn(inpFile);
            k=k+1;
          end 
          
          end % if Flight Data
          
          vInp = 0;
          if strncmpi('Rotor Data',str2find,5) > 0
          % skip comment lines and read three values currently (as of 03/13) 
          % two are dummy values and only last one is used 
           readDmy(inpFile,3);
           itrType     =readLn(inpFile);
          %  read main rotor initial guesses
            readDmy(inpFile,3);
            xh    = readLn(inpFile);
            thtaN = deg2rad(xh);
            A1S = deg2rad(readLn(inpFile));
            B1S = deg2rad(readLn(inpFile));
            xh    = readLn(inpFile);            
            
          % read tail rotor/ rudder values  
            readDmy(inpFile,5);
            lngCntrlTR  =readLn(inpFile);
            vInp     =readLn(inpFile);
          end

          
          readCont= 0;
      end  %     if k > 0
      % end
    end %   if '%'==token(1)
    end % if length(sAux) > 0
    
  end % while
 
else
  fprintf('Error opening file %s %d\',inFileName,inpFile)
            error(' *** PROGRAM TERMINATED ***')
end

if  strncmpi('Flight Data',str2find,6) > 0
  gloVars.rfrnc.xU= rfrncVal(1); gloVars.rfrnc.xW= rfrncVal(2);  gloVars.rfrnc.xQ= rfrncVal(3);
  gloVars.rfrnc.zU= rfrncVal(4); gloVars.rfrnc.zW= rfrncVal(5);  gloVars.rfrnc.zQ= rfrncVal(6);
  gloVars.rfrnc.mU= rfrncVal(7); gloVars.rfrnc.mW= rfrncVal(8);  gloVars.rfrnc.mQ= rfrncVal(9);

  % convert all velocities to m/s
  % vE(1) =  vE(1)/3.6;
  % vE(2) =  vE(2)/3.6;
  % vE(3) =  vE(3)*1.805/3.6;

  
  thtaClmb = -atan(vE(3)/sqrt(vE(1)^2 + vE(2)^2));
  
end


if  strncmpi('Rotor Data',str2find,5) > 0
% if the initial guess for induced velocity is 1.0e6 or  beyond the value of 
% longitudinal tail rotor cyclic is taken to be an elevator deflection 
% that changes the effective angle of incidence of the horizontal stabilizer
% an angle < 0 increases the fuselage attitude angle theta
  if vInp > 0.99e6
     iRigHS=deg2rad(lngCntrlTR);
  end
end
tmpStr = ['End reading ' str2find ' section\n'];
fprintf(tmpStr);
fclose(inpFile);

return 
