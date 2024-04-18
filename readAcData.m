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
% This file reads all aircraft geometry data. This is not a subroutine
% so all data are immediately available in the workspace. Most of the
% data are distributed throughout the whole program via the gloHeader,
% which is a sort of COMMON block (note the all uppercase, all ye 
% FORTRAN dinosaurs..;-)
%
gloHeader

                         %
 readDmy(inpFile,3)      % Fuselage                     
                         %                              
 xh =readLn(inpFile);    % gross  weight                  [lbs ] -?-
 mHeli = xh*lbMass;
 mHeliG= mHeli*gEarth;
 
 fpA=readLn(inpFile)*f2m^2;     % equivalent flat plate area (along x)[ft^2]
 afH=readLn(inpFile)*f2m^2;     % horizontal projected area  (along y)[ft^2] 
 afV=readLn(inpFile)*f2m^2;     % vertical projected area    (along z)[ft^2]
 zcg=readLn(inpFile)*f2m;       % CG height  above waterline  [ft  ]
 xcg=readLn(inpFile)*f2m;       % CG fuselage  station        [ft  ]
 ycg=readLn(inpFile)*f2m;       % CG position rt of butline   [ft  ]
 
 Ixx=readLn(inpFile)*slgft2kgm; % Ixx          [slug ft^2]
 Iyy=readLn(inpFile)*slgft2kgm; % Iyy          [slug ft^2]
 Izz=readLn(inpFile)*slgft2kgm; % Izz          [slug ft^2]
 Ixz=readLn(inpFile)*slgft2kgm; % Ixz          [slug ft^2]
 vfv1=readLn(inpFile);        % downwash  ratio [         ]
 clNfus = readLn(inpFile);    % fuselage lift and..
 cmNfus  = readLn(inpFile);   % moment constant
 
 % if the first token in the next line is a comment no dedicated
 % data file for fuselage coefficents exists and the default in 
 % RotaryWingSimStart is used
 sAux = fgetl(inpFile);
 token=strtok(sAux);
 if '%'==token(1)
   nSkp=2;  % no input, skip next 2 comment lines
   fusDatFileName='';
 else
   [fusDatFlName errMsg]=sprintf('%s',token);
   nSkp=3; % one more input line, skip next 3 comment lines
 end
                             % 
 readDmy(inpFile,nSkp);      % Main rotor
                             % 
 bTl = readLn(inpFile);      % tip loss factor 
 cLftGlo=readLn(inpFile);    % global fudge factor for lift
 cDrgGlo=readLn(inpFile);    % global fudge factor for drag
 nB    =readLn(inpFile);     % number of blades [     ] -?-
 aN    =deg2rad(readLn(inpFile)); % constant coning angle for teetering rotor [degs ]
 rRot  =readLn(inpFile)*f2m; % rotor  radius    [ft   ] -?-
 cHMR  =readLn(inpFile)*f2m; % blade  chord     [ft   ] -?-
 cHtip =readLn(inpFile);     % ratio (blade chord tip)/(blade chord root) [   ]
 xRotT =readLn(inpFile);     % start of tapered section, 1.0 = no taper   [   ]  
 grip  =readLn(inpFile)*f2m;      % blade grip length (root cut out)      [ft ]
 eRot  =readLn(inpFile)*f2m/rRot; % hinge offset            [ft] -?-
 thta1 =deg2rad(readLn(inpFile)); % blade twist             [degs ]-?-
 airfoil=readLn(inpFile);         % blade airfoil No        [     ]
 aLift =readLn(inpFile);          % blade lift curve  slope [     ] -?-
 aLift =aLift*cLftGlo;

 dDrgN  =readLn(inpFile);          % blade drag coefficient  [     ] -?-
 if dDrgN > 0
     dNMR=dDrgN;
 else
     dDrgN=dNMR;
 end
 dDrg1  =readLn(inpFile);          % blade drag coefficient  [     ] -?-
 dDrg2 =readLn(inpFile);          % blade drag coefficient  [     ] -?-

 nRot=0; vTip=0; oM=0;
 [vAux sAux]  =readLn2(inpFile);  % rotational  velocity    [rads/sec] -?-
 if strcmp('rad/sec',sAux) | strcmp('rads/sec',sAux)
    oM   =vAux;
    nRot = oM/2/pi*60;
    vTip = rRot*oM;
 end
 if strcmp('revs/min',sAux) | strcmp('rpm',sAux)
    nRot = vAux;
    oM = 2*pi*nRot/60;
    vTip = rRot*oM;
 end
 if strcmp('fps',sAux)
   vTip = vAaux;
   oM   = vTip/rRot;
   nRot = oM/2/pi*60;
 end
 if ( oM+nRot+vTip < 1.0e-8 )
    fprintf('Input error: omega or nRot or vTip needed\n');
    error(' *** END OF PROGRAM  ***')
 end
 
 nRotMin = nRot;
 oMR = oM*rRot;

 wBlade=readLn(inpFile)*lbMass;   % blade weight            [lbs  ] -?-
 xGBlade=readLn(inpFile);         % blade center of gravity location          
 tipWeight=readLn(inpFile)*lbMass;
 
 iBlade=readLn(inpFile)*slgft2kgm;% flapping moment of  inertia [slug ft^2]
 gamRot= readLn(inpFile);        % rotor Lock number
 hmd  =readLn(inpFile)*f2m;      % hub height  above waterline [ft ]
 lmd  =readLn(inpFile)*f2m;      % hub fuselage  station       [ft ]
 ymd  =readLn(inpFile)*f2m;      % hubposition rt of butline   [ft ]
 iRigMR=deg2rad(readLn(inpFile));% mast incidence, negative forward [deg] -?-
  % if the first token in the next line is a comment no torsional data exist for the rotor
 sAux = fgetl(inpFile);
 [token rest]=strtok(sAux);
 if '%'==token(1)
   twstData=0;
   nSkp=2;  % no input, skip next 2 comment lines
 else
   cGA      = str2double(token)*f2m;           %  distance aerodynamic to center of mass [ft] -?-
   Gmod     = readLn(inpFile)*lbForce*f2m;  % torsional rigidity [lb-ft] -?-
   CmPrfl   = readLn(inpFile);              % blade profile moment coefficient [     ] -?- 
   twstData =1;   
   nSkp=3; % input finished, skip next 3 comment lines
 end
                       % 
 readDmy(inpFile,nSkp) % Tail rotor
                       %
 nBTR=readLn(inpFile);     % number  of blades                 [  ] -?-
 chTR=readLn(inpFile)*f2m; % blade chord                       [ft]
 rRotTR=readLn(inpFile)*f2m; % blade radius                    [ft]
 aLiftTR=readLn(inpFile);    % balde lift curve  slope         [  ] -?-
 dNTR  =readLn(inpFile);     % blade drag coefficient          [  ] -?-
 oMTR=readLn(inpFile);          % rotational  velocity            [rad/sec] -?-
 IbTR=readLn(inpFile)*slgft2kgm;% flapping moment of inertia   [slug ft^2]
 gamRotTR=readLn(inpFile);      % TR Lock number          [  ] 
 delta3TR=readLn(inpFile);      % Delta-3  angle              [deg]
 delta3TR=deg2rad(delta3TR);
 thta1TR=deg2rad(readLn(inpFile));  % blade  twist             [deg]
 trPwrRatio=readLn(inpFile);   % part of of engine power used up by
                               % transmission losses and tail rotor. 
                               %   0.0 < trPwrRatio < 1.0      [  ]
 htd=readLn(inpFile)*f2m;      % hub height  above waterline   [ft]
 ltd=readLn(inpFile)*f2m;      % hub  fuselage  station        [ft]
 ytd=readLn(inpFile)*f2m;      % hub position rt of butline    [ft]
                        % 
 readDmy(inpFile,3)     % Wing
                           % 
aWing=readLn(inpFile)*f2m^2; % area             [ft^2] -?-
bWing=readLn(inpFile)*f2m;   % span             [ft  ] -?-
CLwing=readLn(inpFile);      % CL               [    ]
if CLwing < 0;  CLwngGlo=-CLwing;  else CLwngGlo= 1.0; end
CdMax =readLn(inpFile);      % CdMax            [    ]
if CdMax < 0;  CDwngGlo=-CdMax;   else CDwngGlo = 1.0;  end
cTwing=readLn(inpFile)*f2m;     % tip  chord       [ft  ]
cRwing=readLn(inpFile)*f2m;     % root chord       [ft  ]
eWing=readLn(inpFile);       % wing efficiency factor      [    ]
alpnWing=readLn(inpFile);      % zero lift angle             [deg ]
% alpNw=deg2rad(alpNw);
iRigWing=deg2rad(readLn(inpFile)); % wing rigging angle  [deg ]
% iw=deg2rad(iw);
aliftWing=readLn(inpFile);   % lift curve  slope           [    ]
hwd=readLn(inpFile)*f2m;     % height  above waterline     [ft  ]
lwd=readLn(inpFile)*f2m;     % fuselage  station           [ft  ]
ywd=readLn(inpFile)*f2m;     % position right of butline   [ft  ]
vwv1=readLn(inpFile);        % rotor downwash  ratio       [    ]
detafdalpfw=readLn(inpFile); % fuselage  downwash  ratio   [    ]
                             %
readDmy(inpFile,3)           % Horizontal Stabilizer
                             %
aHoriz=readLn(inpFile)*f2m^2;    % area              [ft^2] -?-
bHoriz=readLn(inpFile)*f2m;      % span              [ft  ]
CLHoriz = readLn(inpFile);       % CL
CDohoriz   =readLn(inpFile);     % CDo               [    ]
aLiftHS =readLn(inpFile);     % lift curve  slope [    ] -?-
alpNh =deg2rad(readLn(inpFile)); % zero lift angle fixed HS  [deg ]
dleOffs = alpNh;                 % offset for hStab geard to longitudinal cyclic
iRigHS=deg2rad(readLn(inpFile)); % angle  of incidence       [deg ]
hhd=readLn(inpFile)*f2m;         % height  above  waterline  [ft  ]
lhd=readLn(inpFile)*f2m;         % fuselage station          [ft  ]
yhd=readLn(inpFile)*f2m;         % position right of butline [ft  ]
qhq=readLn(inpFile);             % Dynamic pressure  ratio   [    ]
vhv1=readLn(inpFile);            % rotor downwash  ratio     [    ]
detafdalpfh=readLn(inpFile);     % fuselage downwash ratio   [    ]
                             %
readDmy(inpFile,3)           % Vertical Fin
                             %
aVert=readLn(inpFile)*f2m^2; % area                        [ft^2] -?-
bVert=readLn(inpFile)*f2m;   % span                        [ft  ]
CLvert=readLn(inpFile);      % maximum Cl                  [    ]
CDovert=readLn(inpFile);     % CDo                         [    ]
aLiftVert=readLn(inpFile);   % lift curve  slope           [    ]
hvd=readLn(inpFile)*f2m;     % height above  waterline     [ft  ]
lvd=readLn(inpFile)*f2m;     % fuselage station            [ft  ]
yvd=readLn(inpFile)*f2m;     % position right of butline   [ft  ]
alpNv=readLn(inpFile);       % zero lift angle             [deg ]
alpNv=deg2rad(alpNv);        %
qvq=readLn(inpFile);         % dynamic pressure  ratio     [    ]
                             % 
readDmy(inpFile,3)           % Control Slopes & Rigging Angles
                             % 
db1mddele   =readLn(inpFile);% long cyclic pitch/inch defl [deg/in  ] 
da1mddela   =readLn(inpFile);% lat cyclic pitch/inch defl  [deg/in  ]
dthetomddelc=readLn(inpFile);% collective pitch/inch defl  [deg/in  ]
dthetotddelp=readLn(inpFile);% tail rotor pitch change/defl[deg/unit]
maxr=readLn(inpFile);        % control until full rudder   [deg/unit]
A1CL=readLn(inpFile);        % rigging angle lateral cyclic [rad]                     
B1CL=readLn(inpFile);        % rigging angle longitudinal cyclic [rad]                
TRCL=readLn(inpFile);        % rigging angle TR collective for zero pedal travel [rad]
                             %
readDmy(inpFile,3)           %  Time Constants
                             % 
tauB =readLn(inpFile);       % tauB, time constant of Bell bar
tauC =readLn(inpFile);       % tauC
K_l  =readLn(inpFile);       % K_l
                             %
readDmy(inpFile,3)           %  Integration Variables
                             % 
nPsi =readLn(inpFile);       % azimut   18 -?-
nStat=readLn(inpFile);       % stations 15 -?-
                             %
readDmy(inpFile,3)           %  Propeller
                             % 
tauProp = deg2rad(readLn(inpFile));  % [deg]
propDia = readLn(inpFile)*f2m;       % [ft ] 
engineRPM =readLn(inpFile);     % engine [rpm] 
engineHP  =readLn(inpFile);     % engine power [hp]
hgtConstPwr = readLn(inpFile);  % height up to which engine power is constant [ft]
hmpr  =readLn(inpFile)*f2m;  % prop height above waterline [ft ]
lmpr  =readLn(inpFile)*f2m;  % prop fuselage station       [ft ]
ympr  =readLn(inpFile)*f2m;  % prop position rt of butline [ft ]
                             %
readDmy(inpFile,3)           %  tip jet data
                             % 
dDrgJet =readLn(inpFile);    % jet drag coefficient
rOutr = readLn(inpFile)*f2m; % [ft]
if rOutr < 0
  rOutr=rRot;
end

rInnr = readLn(inpFile)*f2m; % [ft] 
if rInnr < 0
  rInnr = -rInnr;
  rCntr = 0.5*(rOutr + rInnr);
  dJet  = rOutr-rInnr;
  aJet =  0.25*pi*dJet^2;
else
  rCntr = readLn(inpFile)*f2m; % [ft]
  if (rCntr > 0)
    hJet = readLn(inpFile)*f2m; % [ft]
    aJet = jJet*(rOutr-rInnr);
  else
    dJet= rInnr;
    rInnr = rOutr - rInnr;
    rCntr = 0.5*(rOutr+rInnr);
    aJet= 0.25*pi*(dJet^2);
  end 
end


fclose(inpFile);
