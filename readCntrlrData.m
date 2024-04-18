function [PIDphi, PIDthta, PIDpsi, PIDcllt, CTRL] = readCntrlrData(fileName,str2find, cKG)
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
%         --- Version 1.2 ---
%
% input:
%   fileName    name of controler variables input file
%   str2find    currently not used
%   cKG         constant offset value for controler data
% 
% output:
%   PIDphi      structures containing values for PID
%   PIDthta     controler that operates about the
%   PIDpsi      corresponding axes
%   PIDcllt     PID controler data for collective pitch
%   CTRL        initial values for collective and cyclic
%               pitch as well as throttel setting
%

global pE vE qE aE usephi f2m lbForce
dataDir='./data/';
inFileName= [dataDir fileName '.cnt'];

fprintf('##     Reading ','\n\n%s\n',str2find);
eval(['disp(''# Data   ',inFileName,'   #'')'])

jCount=0;
inpFile = fopen(inFileName,'r');

if ( inpFile < 0 )
  fprintf('\nInput error: no control  file %s \n\n',inFileName);
  error(' *** END OF PROGRAM  ***')
end

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
        fprintf('Start reading control data\n');
        readDmy(inpFile,1)
        c.Kp  = readLn(inpFile) + cKG;
        c.Kd  = readLn(inpFile) + cKG;
        c.Ki  = readLn(inpFile) + cKG;
        readDmy(inpFile,1)
        PIDphi.Kp     = readLn(inpFile)*c.Kp;
        PIDphi.Kd     = readLn(inpFile)*c.Kd;
        PIDphi.Ki     = readLn(inpFile)*c.Ki;
        PIDphi.dlaUpr = readLn(inpFile);
        PIDphi.dlaLwr = readLn(inpFile);
        readDmy(inpFile,1)
        PIDpsi.Kp     =  readLn(inpFile)*c.Kp;
        PIDpsi.Kd     =  readLn(inpFile)*c.Kd;
        PIDpsi.Ki     =  readLn(inpFile)*c.Ki;
        PIDpsi.dlpUpr =  readLn(inpFile);
        PIDpsi.dlpLwr =  readLn(inpFile);
        readDmy(inpFile,1)
        PIDthta.Kp    = readLn(inpFile)*c.Kp;
        PIDthta.Kd    = readLn(inpFile)*c.Kd;
        PIDthta.Ki    = readLn(inpFile)*c.Ki;
        PIDthta.dleUpr= readLn(inpFile);
        PIDthta.dleLwr= readLn(inpFile);
        readDmy(inpFile,1)   
        PIDcllt.Kp    = readLn(inpFile)*c.Kp;
        PIDcllt.Kd    = readLn(inpFile)*c.Kd;
        PIDcllt.Ki    = readLn(inpFile)*c.Ki;
        PIDcllt.zUpr  = readLn(inpFile);
        PIDcllt.zLwr  = readLn(inpFile);
        readDmy(inpFile,1)
        CTRL.dlc  = readLn(inpFile);
        CTRL.dla  = readLn(inpFile);
        CTRL.dle  = readLn(inpFile);
        CTRL.dlp  = readLn(inpFile);

        
        readCont= 0;
      end
    end
    end
  end % while
else
  fprintf('Error opening file %s %d',inFileName,inpFile)
  return
end


fprintf('End read controler data\n');
fclose(inpFile);

return