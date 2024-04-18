%      RotaryWingSimStart
%
% A program for the analysis of rotary wing aircraft dynamics 
%
% Copyright 2010 - 2013 Juergen Humt
% 
%        Version 1.8
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
% The first version of this program was based on naca report naca-tm-73254. 
% The numbers of some of the formulae given in this report are listed in 
% brackets in the code line of the modules that belonged to v 1.0
%
% Revision history       
%                       --- Version 1.0 ---      09/19/10
%                            1.01                10/26/10
%   added main rotor rigging angle (iRigMR) to force components in MainRotorFrcMmtFunc,
%   formula 17 and 19
%
%                            1.2                 05/20/11
%   modified calculation of longitudinal and lateral cyclic in hub wind axes in 
%   MainRotorCtrlFunc. Added negative branch of if statement for small B1S.
%                            06/12/11
%   error in propeller function: prop moment acts about z-axis should be y-axis
%                            06/16/11
%   Qmr was not muliplied by 0.5*sigma, thus rotor accellerations were to high by two 
%   orders of magnitude
%                            07/03/11
%   using drag coefficient from file for vertical fin
%                            08/23/11
%   alfaD = alfaNf + a1, naca-716 uses alfaNf, in some places in the program alfaD 
%   was erroneously used instead
%                            09/17/11
%   wing model added for 1920/30th style gyros. Thrust coefficient is calcu-
%   lated based on force summation, not aircraft mass
%                            12/04/11
%   in the flight data file rpm or oMR can be entered instead of mu. thtaN 
%   is again the first item in the extra section
%                            12/23/11
%   induced velocity used to calculate body velocity for horizontal
%   stabilizer, fuselage and wing. Simulink model thus modified.
%
%                           Version 1.3
%                            15/05/12
%   model cleaned up. File ouput for all contributions to stability derivatives: 
%   rotor, H-stab, V-fin, wing, fuselage. All output data written to file.
%
%                           Version 1.4
%                            28/10/12
%   fuselage lift, drag and moment calculated from polynomial coefficients. 
%   Bramwells curves for calculating induced velocity within the rotor flow field 
%   are used. Induced velocity vi is applied to horizontal stabilizer, wing and 
%   fuselage. Additional data section comprises 7 items #1-> added mass,
%   #6-> thtaN
%                             18/11/12
%   corrected error reading CG offset from flight data file xCG was calculated
%   should have been xcg (case sensitive!)
%
%                           Version 1.5
%                            12/09/12
%   rotor speeds are read from a single line in the aircraft data *.dat file
%   The distinction as to whether the input is revolutions per minute [revs/min], 
%   tip speed in feet per second [fpm] or omega [rad/sec] is by reading a string 
%   enclosde in square brackets which is in the comment following the numerical 
%   value. Before this line there are three lines from which the coefficients for 
%   the three term drag polare are read. The three input lines for blade weight, 
%   blade center of gravity location and tip weight follow the rotor speed 
%   input. Due to the changes the input files are not compatible with any previous 
%   version.
%                           Version 1.6
%                            02/24/13
%   the rotor data section of the flight data file has not been used so far. Now
%   for gyros (no tailrotor) the longitudinal input from the tail rotor is inter-
%   preted as an elevator setting which alters the angle of attack of the hori-
%   zontal stabilizer.
%   the value for lift coefficient of wing is used as a global fudge factor for
%   wing lift if less than zero
%                            02/27/13
%   measured rotor angles a0 a1 b1 and alfaNf are read from file for comparison
%                            03/19/13
%   the program can either search the required collective thtaN at a given rota-
%   tional speed oM or vice versa. The variable controlling the iteration type
%   is read in the first block of data in the *.flt file. The value for thtaN
%   has been moved from the additional data section at the end of the flt file 
%   to the main rotor section at the very beginning. Pedal control input is 
%   interpreted as rudder deflection for autogyros.
%
%                           Version 1.7
%                            04/06/13
%   added engine torque to rotor torque equation. Simulation may be used in 
%   either autogyro or helicopter mode. Induced velocity at the horizontal stabi-
%   lizer vIx is multiplied by rotor downwash ratio vhv1 which is zero if HS is 
%   outside rotor disk. Yaw control is switched between rudder and tail rotor 
%   depending on fligth mode since both may have e.g. different CG distances. 
%   Fuselage moment coefficient cmfNFus is used as a multiplier for the moment 
%   curve read from file to allow quick changes of fuselage contribution.
%                            04/20/13
%   corrected error in calculating dxMR (distance between thrust line and CoG)
%
%                           Version 1.71
%                             04/27/13
%   zero lift angle of horizontal stabilizer is used as offset value for
%   stabilizers which are geared to longitudinal cyclic (dleOffs)
%                             07/31/13
%   horizontal stabilizer gearing curve is read from file if H-stab rigging
%   angle is set to -1.5e-6 (crazy value, eh?...;-)
%   rotor rpm from the *.dat file is overwritten if value in *.flt is > 0
%
%                           Version 1.8
%                             06/10/13
%   added routine for preliminary helicopter performance calculation based on the 
%   paper by Petrovi/Stupar/Kosti/Simonovi, title: "Determination of Light Helicopter 
%   Flight Performance at the Preliminary Design Stage". This calculation is only
%   performed if the variable itrType in the rotor data section of the *.flt file
%   is set to 2
%
%                       --- Version 2.0 ---
%                             11/01/13
%   first release version. Code cleaned (a bit), parameter lists added to all
%   new functions, documentation (a bit) completed.
%
%
%                           Version 2.0.1
%                             23/11/13
% vfv1 added to fuselage force and moment function
%
%                           Version 2.0.2
%                             02/08/15
% rhoAirN introduced to modify rotor Lock number (gamRot) for different
% values of density
%
% General Remarks:
% 
% Capital N in any variable means, that this variable has a zero (German = Null) 
% as subscript. Thus collective pitch is thtaN (theta Null). A capital or lower 
% case v denotes a ratio, thus the varible for P/L is PvL (P versus L)
%
% Several modules had to be added to turn the formulae of the report naca 73254 into an 
% operational mathematical model. Generic 6 DOF body dynamics equations have been added 
% as well as coordinate transformations. For inflow calculation the model uses the ite-
% rative scheme proposed by G. D. Padfield. The trim routine essentially is the one 
% published by M. E. Dreier. Propeller forces and moments were added for rotary wing 
% aircraft such as compound helicopters and autogyros. Wings, as used in many gyros of 
% the 30s, were added as well. Rotor speed was added as 8th degree of freedom  in auto-
% gyro mode and roto torque as 9th degree of freedom in helicopter mode. Fuselage lift,
% drag and moment coefficients are provided as polynomials so measured curves can easily
% be incorporated into the calculation. 
%
% 
%
%  Il semble que la perfection soit atteinte
%   non quand il n'y a plus rien a ajouter,
%  mais quand il n'y a plus rien a retrancher.
%   - Antoine de Saint-Exupèry -
%
%  Perfection is seemingly attained
%  not when there is no longer anything to add,
%  but when there is no longer anything to omit
%
%  Vollkommenheit entsteht scheinbar nicht dann
%  wenn man nichts mehr hinzufuegen kann, sondern
%  wenn man nichts mehr wegnehmen kann.
%
% 
clear all
clc
close all


gloHeader

gloVars.Fx=0; gloVars.Fz=0;

tMax=12;  dt=0.5;  jTOut=1;  jtCoeff= -1; octRun=0;

k=tMax/dt + 1;

nxT=zeros([1  k]);
t=dt; 

for j=1:k
  nxT(j)=t;
  t=t+dt;
end

gEarth  =  9.81;   % acceleration of gravity m/s^2
kinVisc = 18.6e-6; % Kinematic viscosity of air m^2/s

%
% Lower boundary for Reynolds number of a gyro.
% n= 200 rpm -> Vblade=2*pi*n/60= 20.944 [1/s]; c=0.3 [m]; R=6 [m]
% Re = 20.944*6*0.3/18.6e-6 = 2.0e6
%
% ### Unit conversion section. Set usUnits=0 if input is in SI units ###
% This feature is not tested !!
usUnits = 1.0;

if (usUnits > 0.5)                % conversion factors
    f2m       = 0.3048;            % feet to meter 
    i2m       = f2m/12;            % inch to meter 
    slg2kg    = 14.6;              % slugs to kg
    slgft2kgm = 1.3564;            % slug*feet to kg*m
    lbMass    = 0.4536;            % lbMass to kg
    lbForce   = lbMass*gEarth;     % lbForce to N
    rhoAir    = 0.0023762*slg2kg/f2m^3; % slug/ft^3
else
    f2m       = 1.0;
    i2m       = 1.0;
    slg2kg    = 1.0;
    slgft2kgm = 1.0; 
    lbMass    = 1.0;
    lbForce   = 1.0;
    rhoAir    = 1.2252; % slug/ft^3
end

rhoAirN= rhoAir;

epsGlo=1.0e-10; epsGlo3=1000*epsGlo; itMax=150; 

% default directories, do not comment out
airfoilDataDir    ='./data/';  % wing profile data
dataBaseDir1      ='./data/';  % generic data stored here
dataBaseDir2      ='./testCases/';

% replace 1/2 below to switch to desired directory
dataBaseDir = dataBaseDir2;


% the two lines below should not have to be edited
dataDir     = dataBaseDir;
aircraftDataDir      = dataBaseDir;  % aircraft and flight state data


% read measured rotor angles and other values
datFlName ='KD1_mmt';

% ctrlFlName ='cntrlrUH1a1';
ctrlFlName = 'cntrlrS51';

% ~~~~~ Start of Data Section ~~~~~~~
% uncomment desired model or create new one
%
% #-#-#-# OH-6 #-#-#-#
% aircraftDataDir='./data/OH-6/';
% % Routh D 4.835
% rotFileName='OH-6A_1';
% fDatFlName ='fDataOH-6A_6N';
% 
% % Routh D -2.85
% rotFileName='OH-6A_1_60kts_1';
% fDatFlName ='fDataOH-6A_60kts_1';
% fDatFlName ='fDataOH-6A_3N';
% 
% 
%  #-#-#-# UH-1 #-#-#-#
dataDir         =[dataBaseDir 'UH-1/'];
aircraftDataDir =[dataBaseDir 'UH-1/'];
rotFileName='RotDatUH1b';
% 
% fDatFlName ='fDataUH1_3N';
% fDatFlName ='fDataUH1_6N';
fDatFlName ='fDataUH1_1NN2';
% fDatFlName ='fDataUH1_15N';
% fDatFlName ='fDataUH1_1NN_2'; % bestApprx_042713Snap1233';

% 
%  #-#-#-# UH-60 #-#-#-#
% rotFileName='RotDatUH-60a';
% fDatFlName ='fDataUH-60_1NN';
% 

%  #-#-#-# S-51 #-#-#-#
% dataDir         =[dataBaseDir 'S-51/'];
% aircraftDataDir =[dataBaseDir 'S-51/'];
% 
% rotFileName ='RotDatS-51_u';
% fDatFlName  ='fDataS-51_mu20_1p22';
% fDatFlName  ='fDataS-51_mu20d';
% datFlName   ='S-51_angles';
% 

%  #-#-#-# R-6763 #-#-#-#
% dataDir         =[dataBaseDir 'R-6763/'];
% aircraftDataDir =[dataBaseDir 'R-6763/'];
% datFlName   ='R-6763_angles';

% rotFileName ='RotDatR-6763_3';
% fDatFlName  ='fDataR-6763_mu20_3';

% rotFileName ='RotDatR-6763_2_e2';
% fDatFlName  ='fDataR-6763_mu20_e2';
% 
% the file below (fDataR-6763_mu30_n194) gives complete output for 
% all figures with input file RotDatR-6763_2_e2
% fDatFlName  ='fDataR-6763_mu10_n194';
% fDatFlName  ='fDataS-51_mu20d';
% 
%  #-#-#-# S-58 #-#-#-#
% aircraftDataDir     ='./data/S-58/';
% rotFileName ='RotDatS-58';
% fDatFlName  ='fDataS-58_mu20';
% datFlName   ='S-51_angles';
% 

%  #-#-#-# R4B #-#-#-#
% aircraftDataDir  =[dataBaseDir 'R4-B/'];
% rotFileName= 'RotDatR4B_Fa';
% % rotFileName= 'RotDatR4B_1a';
% 
% 
% jF=  6;
% if jF < 10
%   xHS=sprintf('%1d',jF);
% else
%   xHS=sprintf('%2d',jF);
% end
% fDatFlName =strcat('fDataR4B_no',xHS);
% 
% 
% fDatFlName ='fDataR4B_40mph';
% fDatFlName ='fDataR4B_no12_a1_12p3';
% fDatFlName ='fDataR4B_c1NN';


%  #-#-#-# R-22 #-#-#-#
% rotFileName='RotDatR22t';
% fDatFlName ='fDataR22_60mph';


% #-#-#-# PCA-2 #-#-#-#
% naca 434
% dataDir         =[dataBaseDir 'PCA/'];
% aircraftDataDir =[dataBaseDir 'PCA/'];
% 
% rotFileName='PCA2tws';
% fDatFlName ='fDataPCA-2_A12_4'; % ++
% fDatFlName ='fDataPCA-2_lvl'; % --
% fDatFlName ='fDataPCA-2_A12_7'; % ++
% fDatFlName ='fDataPCA-2_A5_5'; % ++
% fDatFlName ='fDataPCA-2_A5_8'; % ++
% fDatFlName ='fDataPCA-2_A6_3'; % ++
%
% naca 475
% fDatFlName ='fDataPCA-2_475_5'; % ++
% fDatFlName ='fDataPCA-2_475_10';% ++
% fDatFlName ='fDataPCA-2_475_18';% ++


% #-#-#-# Prouty #-#-#-#
% iTyp=1;
% if (2==iTyp)
%   rotFileName='RotDatProuty2';
%   fDatFlName ='fDataPrty-2';
% else
%   rotFileName='RotDatProuty6';
%   fDatFlName ='fDataPrty-7';
% end
% 

%  #-#-#-# Puma  #-#-#-#
% dataDir         =[dataBaseDir 'Puma/'];
% aircraftDataDir =[dataBaseDir 'Puma/'];
% rotFileName='RotDatPuma1b';
% rotFileName='RotDatPuma1a_phgd_p_shrtTrm_Ok';
% fDatFlName ='fDataPuma14Nkt';
% fDatFlName ='fDataPuma4Nkt';
% rotFileName='RotDatPuma1c';
% fDatFlName ='fDataPuma1c';
% 

% 
% #-#-#-# XV-1 #-#-#-#
% rotFileName='XV-1a';
% fDatFlName ='fDataXV-1_1';

% #-#-#-# Montgomerie Parsons #-#-#-#
% rotFileName = 'MntGmPrZ';
% rotFileName = 'MntGmPrNoCm';
% rotFileName = 'MntGmPrNoCm_smallB1C';

% fDatFlName  = 'fDataMntGmPr_35';
% fDatFlName  = 'fDataMntGmPr_4N';
% fDatFlName  = 'fDataMntGmPrNoCm';
% fDatFlName  = 'fDataMntGmPrNoCm_smallB1C';


% #-#-#-# Bensen #-#-#-#
% dataDir         =[dataBaseDir 'Bensen/'];
% aircraftDataDir =[dataBaseDir 'Bensen/'];
% 
% rotFileName='RotDatBensenB8M_a';
% fDatFlName ='fDataBensenB8M_a';
% fDatFlName ='fDataBensenB8M_60mph';
%  

% #-#-#-# VPP M-16 #-#-#-#
% dataDir         =[dataBaseDir 'VPP-M16/'];
% aircraftDataDir =[dataBaseDir 'VPP-M16/'];
% 
% rotFileName='VPMM16t';
% fDatFlName ='fDataVPMM16a2';
% 


% #-#-#-# Kellet KD-1 #-#-#-#
% Kell..tws  -> rotor speed de(!)creases with flight speed
%
% Kell..tws3 -> rotor speed in(!)creases with flight speed
% rotFileName= 'KelletKD1wtws5';
% rotFileName= 'KelletKD1tws3';
% fDatFlName = 'fDataKD1tws4';
% 

%  #-#-#-#  HS-1x  #-#-#-#
% aircraftDataDir     ='./data/HS-1x/';
% rotFileName='RotDatHS-1x_3';
% fDatFlName ='fDataHS-1x_3';

% ~~~~~ End of Data Section ~~~~~~~



logFileName= [dataDir rotFileName '.log'];
logFile = fopen(logFileName,'w');


outFileName= [aircraftDataDir fDatFlName '.lis'];
outFile = fopen(outFileName,'w');


inFileName= [aircraftDataDir rotFileName '.dat'];
inpFile = fopen(inFileName,'r');
fprintf(outFile,'reading aircraft data from %s\n',inFileName);
fprintf('reading aircraft data from %s\n',inFileName);

% The file readAcData is an include rather than being a function. It takes
% no parameters (i.e it is no subroutine) and therefore all data read in 
% are immediately available in the workspace.
readAcData;

zB  = rRot/1.35;
G1  = 0.2;
kG  = 1-exp(-zB/(2*rRot)/G1); % ground effect factor
lH = (1 - 0.5/(1 + 4.0*(zB/rRot)^2));


% if rotor twist data were read the formulae of
% naca 600 are used for rotor calculation
% twstData = 0;
if twstData < 0.5
  clcCT_hndl     = @clcCT;
  clcRtAngl_hndl = @clcRotAngl;
  fprintf('Twist NOT included\n');

else
  clcCT_hndl     = @clcCTtwst;
  clcRtAngl_hndl = @clcRotAnglTwstKx;  
  fprintf('Twist IS included\n');
end


% values for three term drag polar cDrg[N,1,2] are
% either read from file or default values are set
if dDrgN < 1.0e-5
  [dDrgN dDrg1 dDrg2] = dDrg;
end

% if no file name for fuselage data was read from the
% aircraft data file a default is used
if length(fusDatFlName) < 1
  fusDatFlName ='PrtyCoeff';
end

% read fuselage lift drag and moment coefficients
inFileName= [aircraftDataDir fusDatFlName '.cff'];
inpFile = fopen(inFileName,'r');

if inpFile > 0
  eval(['disp(''# fuselage coefficients   ',inFileName,'   #'')'])
  [ClFus] = cXstrct(inpFile);
  [CyFus] = cXstrct(inpFile);
  [CdFus] = cXstrct(inpFile);
  [CmFus] = cXstrct(inpFile);
  [CnFus] = cXstrct(inpFile);
else
  fprintf('Error opening coefficient file %s.cff\n using generic coefficients\n',fusDatFlName)
  return
end
fclose(inpFile);

CmFus.cAmp = CmFus.cAmp + cmNfus;

inFileName= [aircraftDataDir datFlName '.cff'];
inpFile = fopen(inFileName,'r');

if inpFile > 0
  eval(['disp(''# rotor angle coefficients   ',inFileName,'   #'')'])
  [aNCff] = cXstrct(inpFile);
  [a1Cff] = cXstrct(inpFile);
  [b1Cff] = cXstrct(inpFile);
  [laCff] = cXstrct(inpFile);
  [cTCff] = cXstrct(inpFile);
  [rrpmCff] = cXstrct(inpFile);    
else
  fprintf('Error opening coefficient file %s.cff\n using generic coefficients\n',datFlName)
  return
end
fclose(inpFile);

% mu=0.225;
% xh = crvVal(mu, laCff);
% 

% read cL vs thtaN[rad] and cD vs cL for power calculation
datFlName = 'cL_alf';
inFileName= [dataDir datFlName '.cff'];
inpFile = fopen(inFileName,'r');

if inpFile > 0
  eval(['disp(''# cL cD coefficients for power  ',inFileName,'   #'')'])
  [cLpwrCff] = cXstrct(inpFile);
  [cDpwrCff] = cXstrct(inpFile);
else
  fprintf('Error opening coefficient file %s.cff\n using generic coefficients\n',datFlName)
  return
end
fclose(inpFile);
 
% read lift drag and moment data for the wing
airFoilStr='NACAM3MR1pNe5';
inFileName= [airfoilDataDir airFoilStr '.pol'];
[eF] = readAirfoilData(inFileName,2);
if ( 1 == eF)  return;  end

% set values for propeller thrust polynomial
pCTP = setThrustPoly;

[thtaNfd alfaD mu]=readFlightData(aircraftDataDir,fDatFlName,'Flight Data',outFile);
[thtaNfd alfaD mu]=readFlightData(aircraftDataDir,fDatFlName,'Rotor Data',outFile);
%fprintf(outFile,'reading flight data from %s\n',fDatFlName);

vFlight=norm(vE);

lWN     = lwd - xcg; % CG offset
hWN     = hwd - zcg; % 

pB = 0; qB = 0; rB = 0.0;
omB = [pB qB rB];

vF= norm(vE);
fprintf(outFile,'\n\n    #### RotaryWingSim ####\n');
fprintf(outFile,'rhoAir       %7.4f kg/m³\n',rhoAir);

zetaSet(0);
dleSet(0);


% ### Main Rotor parameter section ###
fprintf(outFile,'thtaN %7.2f°   twist  %7.2f°\n',rad2deg(thtaN),rad2deg(thta1));
xCG = lmd - xcg; % center of gravity (CofG) offset to main rotor
hMR = hmd - zcg; % 2.07/;  naca-tm-73254 uses lh


aDiskMR = pi*rRot^2; % rotor disk Area
diskLoad= mHeli/aDiskMR;
% fprintf(outFile,'dskLd        %7.2f kg/m^2\n',diskLoad);
fprintf(outFile,'mHeli          %7.2f lb\n',mHeli/lbMass);

if 2==itrType
fprintf(outFile,'\nnRot   %8.1f 1/min    oM   %7.2f rad/s\n',nRotMin,oM);
fprintf(outFile,'oMR     %8.2f m/s     oMR %8.2f ft/s\n',oMR,oMR/f2m);
fprintf('\nnRot   %8.1f 1/min    oM   %7.2f rad/s\n',nRotMin,oM);
fprintf('oMR     %8.2f m/s     oMR %8.2f ft/s\n',oMR,oMR/f2m);
end


sigma = nB*cHMR/pi/rRot;          % rotor solidity
% bTl= 0.970;   xNr=0.0;          % tipp loss factor and root cutout
cTCx = 0.5*aLift*sigma;
% aDisk  = pi*rRot^2;
aBlade = aDiskMR*sigma;
fprintf(outFile,'aBlade     %9.2f m^2\n',aBlade);


% ### tail rotor Parameter Section ###
lTR = ltd - xcg; % 8.53;  distance from C of G to TR
hTR = htd - zcg; % 2.1;  

aDiskTR= pi*rRotTR^2;

oMRTR= oMR*rRotTR; % ft/s
n=oMTR/(2*pi);
nRotMinTR=n*60;    % 1/s rotational speed (rpm)


sigmaTR  = nBTR*chTR/pi/rRotTR;          % tail rotor solidity

% output of tailrotor data is suppressed if 
% radius of tailrotor is less than 0.1 [m]
outTr=10.0*rRotTR;
if outTr > 0.5 
fprintf(outFile,'\nrRotTR %10.2f m\n',rRotTR);
fprintf(outFile,'%s%10.1f%s\n','nRotTR   ',nRotMinTR,' 1/min');
fprintf(outFile,'oMTR  %8.2f\n',oMTR);
fprintf(outFile,'oMRTR %8.2f m/s    oMRTR %8.2f ft/s\n',oMR,oMR/f2m);
fprintf(outFile,'sigmaTR %8.4f\n\n',sigmaTR);
end

if (usUnits < 0.5)
  fprintf(outFile,'%s%9.2f%s\n','Area Mrot  ',aDiskMR,' ft^2');
else
  fprintf(outFile,'%s%9.2f%s\n','Area Mrot  ',aDiskMR,' m^2');
end


hFlight    = 0.0;

if 2==itrType
  helicopterPerformance (hFlight, oM, trPwrRatio, hgtConstPwr);
end

% return
% ========== heli Parameter Section ==============

% ### Fuselage Parameter Section ###
% Equivalent flat plate area scaled by ratio of rotor radius.
% This line may be used if no fpA value is available and gives
% an estimate for a UH1-H type helicopter. For a fairly draggy 
% helicopter (e.g. R4-b) replace the 1.0 by about 2.0
% fpA2   =  1.0*(rRot/2.125);
% 
mHeli  = mHeliG/gEarth;

InrtMatHeli = [Ixx 0 Ixz; 0 Iyy 0; Ixz 0 Izz];

Sref  = 1.32*fpA; % (48.0) these constants are calculated..
lref  =  1.4*lTR;  % (39.0) ..backwards for the UH1-H


% ### Blade Parameter Section ###
% blade parameters are approximately calculated if they are not provided on 
% file using the mass constant of the rotor blade gamRot= c*ro*a*r^4/I1
%
roAL  = 2700;          % kg/m1^3 density of aliminium alloy
tRot  = 0.000450*rRot;  % mm thickness of blade skin based on rRot  


if (wBlade < 1.0e-8)
  wBlade = 2*tRot*cHMR*rRot*roAL;     % kg mass of rotor blade
  fprintf(outFile,'%s%10.1f%s\n','m blade  ',wBlade,' kg     estimated');
else
  fprintf(outFile,'%s%10.1f%s\n','m blade  ',wBlade,' kg');
end  

if (iBlade < 1.0e-8)
  iBlade = 1/3*wBlade*rRot^2;      % kg*m^2 moment of inertia about hinge;
  fprintf(outFile,'I blade  %10.1f kg*m^2 estimated\n',iBlade);
else
  fprintf(outFile,'I blade  %10.1f kg*m^2\n',iBlade);
end  

if (gamRot < 1.0e-8)
  gamRot  = cHMR*rhoAir*aLift*rRot^4/iBlade;
end


mWRot = wBlade*rRot/2*gEarth;   % weight moment about hinge
cMIN  = mWRot/iBlade; % constant used in naca-tr-716

% gamRot has to be corrected for specific density of air
gamRot = gamRot*rhoAir/rhoAirN;
fprintf(outFile,'gamRot      %10.4f\n',gamRot);
fprintf(outFile,'sigma       %10.4f\n\n\n',sigma);

% fprintf(outFile,'thtaN %6.2f\n',rad2deg(thtaN));

% left and right limit for mu search
xL=0.001; xR=0.75;

% uncomment the line "clcNACA2154" below to call the program segment 
% that calculates glide performance using the formulae of naca 2154
% clcNACA2154


% ### ------ Flight Angles Trim Values -------- ###
% thtaTrim is a first guess 
phiTrim= aE(1);  thtaTrim= aE(2);  psiTrim= aE(3);

fprintf(outFile,'Trim phi %10.2f°    thta %10.2f°    psi %10.2f°   \n',rad2deg(phiTrim),rad2deg(thtaTrim),rad2deg(psiTrim));
fprintf(outFile,'vH %10.2f mph     vV %10.0f fpm\n',vE(1)*3.6/1.605,vE(3)*60/f2m);
xCTNN = rhoAir*aDiskMR*(oM*rRot)^2;
cTReq = mHeliG*cos(thtaTrim)/xCTNN; %/cos(deg2rad(3.7));
cTsN  = cTReq/sigma;

muInp=norm(vE)/oMR*cos(deg2rad(3.5));  alfaNf=deg2rad(5);
cL=2*cTReq/muInp^2*cos(alfaNf)^3;
fprintf(outFile,'cTs %9.6f     cT %7.4f     cT2sa %7.4f     cL %7.4f\n',cTsN,cTReq,2*cTReq/sigma/aLift,cL);
fprintf('cTs %12.8f     cT %7.4f     cT2sa %7.4f     cL %7.4f\n',cTsN,cTReq,2*cTReq/sigma/aLift,cL);
fprintf('rRot %8.2f     sigma %7.4f\n',rRot,sigma);

% ##### ======== ########
%

% Constants used in report naca-cr-73254
  R3=0.5/aLift;  % R3=0.087;
  R5=0.25/aLift; % R5=4.36e-2;

  kB=K_l*tauB; % 0.528

  C1 =db1mddele   ; % Long. main rotor cyclic angle vs stick TRavel [rad/cm]
  C4 =da1mddela   ; % Latr. main rotor cyclic [rad/cm]
  C5 =dthetomddelc; % main rotor collective [rad/cm]
  C6 =dthetotddelp; % TR collective vs pedal travel [rad]
  C7 =TRCL;         % rigging angle TR collective for zero pedal travel [rad]

  RNTR = rhoAir*pi*rRotTR*oMRTR^2;
  T1 = 0.5/oMRTR;   %  22.17e-4
  T2 = sqrt(RNTR)*aLiftTR*sigmaTR/(8*sqrt(2)); % 35.7
  T3 = RNTR*aLiftTR*sigmaTR/6.51;  %  31.10e+3 /6 was adjusted to /6.51
  T4 = 217.4; % *rRotTR/1.5;
  T5 = 0;
  
  % fuselage drag constants from projected areas
  D1 = 0.5*rhoAir*fpA;
  D2 = 0.5*rhoAir*afH;
  D3 = 0.5*rhoAir*afV;
  
  L1 = 0.5*rhoAir*Sref*clNfus;   % 0.319
  Y1 = L1;

  M1a = 0.5*rhoAir*Sref*lref*cmNfus; % 13.95
  M1 = 0.5*rhoAir*rRot^3;
  N1 = L1;

if (usUnits < 0.5)
  fprintf(outFile,'m Heli %10.1f lb     fpA %10.1f ft^2\n',mHeliG,fpA);
  fprintf(outFile,'%s%10.1f%s\n','Iyy Heli    ',Iyy,' slugs*ft^2');
else
  fprintf(outFile,'m Heli %10.1f N      fpA %10.1f m^2    %10.1f ft^2\n',mHeliG,fpA,fpA/f2m^2);
  fprintf(outFile,'%s%10.1f%s\n','Iyy Heli    ',Iyy,' kg*m^2');
end

% ### ---- wing section ----- ###
lWN     = lwd - xcg; % CG offset
hWN     = hwd - zcg; % 

% ### ---- horizontal stabilizer section ----- ###
areaHS  = aHoriz; % 16.4*f2m^2; % HS area
bHS     = bHoriz; % 8.75*f2m;   % span

lHS     = lhd - xcg; % 6.0  m
hHS     = hhd - zcg; % 0.5  m;

vHS     = lHS*rRot*areaHS/(sigma*aDiskMR);
% fprintf(outFile,'%s%10.4f\n','vHS      ',vHS);

H1 = 0.5*rhoAir*areaHS*aLiftHS;   % 1.48
H2 = H1*tan(deg2rad(5));

CDalf9N=0.5;
H4 = 0.5*rhoAir*areaHS*CDalf9N;   % 0.4

% ### |||| Vertical Stabilizer Section |||| ###
areaVSfin  = aVert; % vertical stabilizer (or fin) area
bVS        = bVert; % 1.37 vertical fin span
aLiftVSfin = aLiftVert; % 5.500 lift curve slope of vertical fin
lVF = lvd - xcg; % 7.62;
hVF = hvd - zcg; % 

vVS = lVF*rRot*areaVSfin/(sigma*aDiskMR);

F1 = 0.5*rhoAir*areaVSfin*CDalf9N;  % 0.77
k1 = 0.5*rhoAir*areaVSfin*aLiftVSfin;    % 1.4
k2 = k1*tan(deg2rad(20));      % 0.51
% fprintf(outFile,'%s%10.4f\n','vVS      ',vVS);


% ### engine/propeller section ###
lPR=lmpr-xcg;
hPR=hmpr-zcg;

fprintf(outFile,'Coordintats x,z    [m]         [ft]\n');
coorOut('Rotor      ',xCG,hMR,outFile);
coorOut('Rudder     ',lTR,hTR,outFile);
coorOut('Fin        ',lVF,hVF,outFile);
coorOut('H-Stab     ',lHS,hHS,outFile);
coorOut('Wing       ',lWN,hWN,outFile);
coorOut('Propeller  ',lPR,hPR,outFile);


% dlpSlct is used to select either tail rotor or rudder for yaw control
% dlpSlct=1-> rudder  dlpSlct=2-> tail rotor
if itrType < 2
  TRCL= 0; % tail rotor rigging anle zero
  C7= 0;
  C8= 0;
  T2= 0;
  dlpSlct= 1;
else
  C8= 0.0873;
  dlpSlct= 2;  
end
dlpArry=[0 0];

% vB = [30 0 -3];
% [fPR, mPR] = PropellerFrcMmtFunc( vB, xh);
%
% ### flight variable section ###
fprintf(outFile,'%s\n','### Var ###');

fprintf(outFile,'mu         %10.3f  vFlight %8.1f kn\n',mu,vFlight*3.6/1.805);

% muTR = vFlight/oMRTR;
% alfaDTR=0.02;
% thtaNTR=0.15;
% cTsTR  =0.068;

% time constant to inhibit algebraic loop
tauAgLp=1.0e-12;


% ### ------ Trim Section -------- ###
% ground effect is currently not used
zB= 10000;

pB = 0; qB = 0; rB = 0.0;

% vErth includes climb/descend anlge 
omDotB = [pB qB rB];
vOut = RTrfFunc(phiTrim, thtaTrim,  psiTrim, vE(1), vE(2), vE(3));

uWStrt = vOut(1); vWStrt = vOut(2); wWStrt = vOut(3);
xDotW = [uWStrt vWStrt wWStrt];


% dlt = 1.5;
% tProp =thrustCalcUS(vFlight, dlt);

% read controler data and initial settings of aircraft controls
% use cKp to switch on and off the pid controlers. 
cKp= 0.0;
[PIDphi, PIDthta, PIDpsi, PIDcllt, CTRL] = readCntrlrData(ctrlFlName,'Controler Data',cKp);
% PIDpsi.Kp=1.0e-2;  PIDpsi.Ki=2.5e-4;

fprintf(outFile,'%s%10.4f%s%10.4f%s\n','thtaClmb   ',thtaClmb,'  ',rad2deg(thtaClmb),'°');

% control movement
%   main         and   tail collective
%   0 == 27.9         17.5L    0.0  17.5R
%                      18° = +4.0 = -10°
dlc= CTRL.dlc;   dlp= CTRL.dlp;
% longitudinal     / lateral cyclic
% -16.38 == 16.38     -16.0 == 16.0
dle= CTRL.dle;   dla= CTRL.dla;
dlp=-0.7;

dlc=thtaN/C5;
xh = rad2deg(thtaN);
% thtaNx=dlc*C5;
% thtaTrim assumes a1s=0
% thtaTrim= thtaClmb -iRigMR + alfaD;



% calculate longitudinal control B1
% B1S = -alfaD + thtaTrim -iRigMR + a1MR;
oMTrim = oM;

% inverse calculate dle from B1S
B1STrim = B1S;
dle=(B1S-B1CL)/C1;

AB1CP  =d2ABCP(dle,dla);

% control input angle
% to give the control plane (swashplate) angles
A1S= AB1CP(1) + A1CL;  A1STrim = A1S;  % AB1CP(1)=   - 2.8
B1t= AB1CP(2) + B1CL;

% cos(alfaD);
cTsReq=cTReq/sigma;       cT2saReq=2*cTsReq/aLift;
fprintf(outFile,'cTsReq %10.5f    cTReq %10.5f    cT2saReq %10.5f\n',cTsReq,cTReq, cT2saReq);

outP=2.0; firstOutC=1;  firstOutB=-1; dlt = 0.0;

% psiTrim=0;
fprintf(outFile,'alfaD %10.2f°    thtaClmb %10.2f°    thtaTrim %10.2f\n',rad2deg(alfaD),rad2deg(thtaClmb),rad2deg(thtaTrim));
% [fR, mR, oMDotMR, fZMR] = AcFrcMmtSummationOutF( dlc, dle, dla, dlp, dlt, phiTrim, thtaTrim, oMTrim, vE, qE, mHeliG, outP);


% return
% set initial values for trim routine
trmVr(1)=dlc;
trmVr(2)=dle;
trmVr(3)=dla;
trmVr(4)=dlp;
trmVr(5)=dlt;
trmVr(6)=phiTrim;
trmVr(7)=thtaTrim;
% trmVr(8)=psiTrim;
trmVr(8)=oMTrim;
trmVr(9)=qMR;

% search trimmed flight values 
jOut   = 1.0;  cDmp= 0.72;
% vE_ar=[40:10:70];
% 
% for jAr = 1:length(vE_ar)
%     vE(1)=vE_ar(jAr);

trmOut = acTrim3DFM(trmVr, vE, omDotB, jOut, cDmp, outFile);
if norm(trmOut) < 1.0e-10
  return
end
% end

% clcDerivs(trmVr, vErth, omDotB,jOut);

dlc= trmOut(1);
dle= trmOut(2);
dla= trmOut(3);
dlp= trmOut(4);  %  -1.3641e-3;
dlt= trmOut(5);
phiTrim  =  trmOut(6);
thtaTrim =  trmOut(7);
oMTrim   =  trmOut(8);
fZMRTrim = -trmOut(9);


jh= -1; outP=-1;
[fR, mR, oMDotMR, fZMR] =   [fext1 mext1 QR1]= AcFrcMmtSummationOutF(trmOut(1), trmOut(2), trmOut(3),trmOut(4), trmOut(5), trmOut(6),  trmOut(7), trmOut(8), trmVr(9), vE, qE, trmOut(9), outP, jh);


if itrType < 2
jOut = 3;
vE2 = vE; vE2(1)=vE2(1)+epsGlo;
trmOut = acTrim3DFM(trmVr, vE2, omDotB, jOut, cDmp, outFile);
DoMDu= (trmOut(8)-oMTrim)/epsGlo;


vE2 = vE; vE2(3)=vE2(3)+epsGlo;
trmOut = acTrim3DFM(trmVr, vE2, omDotB, jOut, cDmp, outFile);
DoMDw= (trmOut(8)-oMTrim)/epsGlo;

omDotB(2)= epsGlo;
trmOut = acTrim3DFM(trmVr, vE, omDotB, jOut, cDmp, outFile);
DoMDq= (trmOut(8)-oMTrim)/epsGlo;

fprintf('DoMDu  %10.6f    %10.6f\n',DoMDu,gloVars.DoMDu_it);
fprintf('DoMDw  %10.6f    %10.6f\n',DoMDw,gloVars.DoMDw_it);
fprintf('DoMDq  %10.6f    %10.6f\n',DoMDq,gloVars.DoMDq_it);
fprintf('DoM_it %10.6f    %10.6f\n',gloVars.DoM_it); 


xUm= gloVars.DxDoM*gloVars.DoMDu_it;    xWm= gloVars.DxDoM*gloVars.DoMDw_it;  xQm= gloVars.DxDoM*gloVars.DoMDq_it;
zUm= gloVars.DzDoM*gloVars.DoMDu_it;    zWm= gloVars.DzDoM*gloVars.DoMDw_it;  zQm= gloVars.DzDoM*gloVars.DoMDq_it;
mUm= gloVars.DmDoM*gloVars.DoMDu_it;    mWm= gloVars.DmDoM*gloVars.DoMDw_it;  mQm= gloVars.DmDoM*gloVars.DoMDq_it;


    fprintf('\n      oM deriv\n');
    fprintf('       u                w              q\n');
    fprintf('x  %10.6f    %10.6f   %10.6f\n',xUm,xWm,xQm);
    fprintf('z  %10.6f    %10.6f   %10.6f\n',zUm,zWm,zQm);
    fprintf('m  %10.6f    %10.6f   %10.6f\n',mUm,mWm,mQm);
    
end


dlcTrim = dlc;
thtaN   = dlc*C5;

AB1CP  =d2ABCP(dle,dla);


% control input angles and rigging angles are added
% to give the control plane (swashplate) angles
A1STrim= AB1CP(1) + A1CL;
B1STrim= AB1CP(2) + B1CL;

fprintf(outFile,'dlc   %6.2f   dla  %6.2f   dle  %6.2f      dlp  %6.2f      dlt  %6.2f\n',dlc,dla,dle,dlp,dlt);
fprintf(outFile,'thtaN %6.2f°   A1S %6.2f°   B1S %6.2f°   thtaTR %6.2f°\n',rad2deg(thtaN),rad2deg(A1STrim),rad2deg(B1STrim),rad2deg(C6*dlp+C7));
fprintf(outFile,'TrimAC phi LFr %10.2f°    thta %10.2f°    psi %10.2f°\n',rad2deg(phiTrim),rad2deg(thtaTrim ),rad2deg(psiTrim));

fclose(outFile);
% ### test section ###
xh = 80*f2m;
vE = [xh,0,0];

lamNf = gloVars.lamNf;
mu    = gloVars.mu;


coeff716(mu);

DvLp = 0.5*rhoAir*norm(vE)^2*fpA/mHeliG;
DvLc = tan(thtaClmb);

cT2sa = eq6(mu, lamNf, thtaN);
cQ2s_1 = (eq11d(mu, lamNf, thtaN) - eq9a(mu, lamNf, thtaN))/aLift;

PvL_1  = cQ2s_1/(mu*cT2sa);



% 0.086
DvLN = abs(eq13(mu,lamNf,thtaN)/(mu*cT2sa));

cTs = 0.5*aLift*cT2sa;
cT  = cTs*sigma;

cNN = (2*mu*(mu^2+lamNf^2)^0.5);
DvLi = cT/cNN;

tanAlfaNf = atan(lamNf/mu+DvLi);

alNf = atan(tanAlfaNf);

alNfDeg = rad2deg(alNf);

cL   = 2*cT*cos(alNf)^3/mu^2;
DvLi_2 = 0.25*cL;

DvLr = DvLN + DvLi;
cLs = cL/sigma;

PvL = DvLN + DvLi + DvLc + DvLp;

DvLpX= DvLp*mHeliG;
xh9  = DvLp*mHeliG;

DvLr = (DvLN + DvLi)*mHeliG;
xh19  = 1.5*DvLN*mHeliG;


xh = PvL/PvL_1;
xh2= PvL_1/PvL;

xh3  = PvL*mHeliG*norm(vE)/1000;
xh3a = xh3/0.734;


fprintf('DvLr  %9.5f = DvLN %9.5f + DvLi %9.5f\n',DvLN+DvLi,DvLN,DvLi);
fprintf('mu    %9.5f   cTs  %9.5f   cT   %9.5f\n',mu,cTs,cT);
fprintf('tht75 %9.2f   cLs  %9.5f   cL   %9.5f\n',rad2deg(thtaN + 0.75*thta1),cLs,cL);

return
