function helicopterPerformance (hFlight, oM, trPwr, hgtConstPwr)
% 
% Input:
% hFlight: height above sea level, meter
% oM     : main rotor rotational velocity
% trPwr  : fraction of engine power used up by 
%          tailrotor and transmission
% tail rotor power as fraction of available power
% hgtConstPwr: height up to engine power is 
%              nearly constant
% 
% capital H in any variable indicates value at height > 0
%

global f2m


% outVal = pwrCalc(0, mu, oM, trPwr, 0, 0);
% vClmb= outVal(1); 
% vDsc = outVal(2);
% vHn  = outVal(3);

% the mu range has to be large enough to cover the maximum 
% flight speed, which often is not at sea level.
muX= [0.0:0.0005:0.35];
vClmb=[];  vHn=[]; vDsc=[];  hgt=[]; cTH=[]; pReq=[];
kC1= 0;  kC2=0; kC2Start=0;
mIndx=[]; mIndx2=[];

vMaxN    = -1.0e5;
vMaxH    = -1.0e5;
clmbStrt =  0.0;

hgtMax   = -1.0e5;
hgtMaxHvr   =0;
hgtMaxHvrIG =0;

h1=0;  h2= 6500;


options = optimset('MaxFunEvals',500,'TolFun',1.0e-2);
% Matlab
hgtHvrIG = fMinbnd(@pwrCalc,h1,h2,options,0, oM, trPwr, hgtConstPwr, 1.2);
% Octave
% hgtHvrIG = fminbnd(@(h1) pwrCalc(h1, 0, oM, trPwr, hgtConstPwr, 1.2) ,h1,h2,options);

% main loop, clalculates maximum rates of climb and descent and uses
% interval search (fminbnd) to determine maximum height for given ad-
% vance ratio (mu= rotorSpeed/flightSpeed)
for j=1:length(muX);
   outVal = pwrCalc(hFlight, muX(j), oM, trPwr, 0, 0);
   vClmb(j) = outVal(1); 
   vDsc(j)  = outVal(2);
   vHn(j)   = outVal(3);
   pReq(j)  = outVal(5);

   % Matlab
   hgt(j) = fMinbnd(@pwrCalc,h1,h2,options,muX(j), oM, trPwr, hgtConstPwr, 0.5);
   % Octave
   % hgt(j) = fminbnd(@(h1) pwrCalc(h1, muX(j), oM, trPwr, hgtConstPwr, 0.5) ,h1,h2,options);
   outVal = pwrCalc(hgt(j), muX(j), oM, trPwr, 0, 0);
   vClmbH = outVal(1); 
   vTmp   = outVal(3);
   cTH(j) = outVal(4);
   
   % the kC1 + kC2 condition eliminates values outside the performance envelope
   if ((abs(hFlight) < 1.0e-5) & (vHn(j) > vMaxN) & ((kC1 + kC2)< 0.5))
      vMaxN = vHn(j);
   end
   
   h1 = hgt(j) - 100;
   if h1 < 0 
      h1 = 0;
   end
       
   h2 = hgt(j) + 100;
   % recording starts after vClmb has been greater zero at least once
   if (vClmb(j) > 0)
     clmbStrt=1;
   end
   
   if (kC1 > 0.5) & (hgt(j) < hgtConstPwr + 1)
      kC2Start =1;
   end
   
   if kC2Start > 0.5
           kC2 = kC2 + 1;
      mIndx2(kC2) = j;
   end
   
   if ((kC1 + kC2)< 0.5) & (vTmp > vMaxH)
       vMaxH = vTmp;
       hMaxV= hgt(j);
   end
   
   if (vClmbH < 0) & (kC2Start < 0.5) & (vTmp > 20)
     if hgtMax  < hgt(j)
       hgtMax = hgt(j);
     end
   end
   
   if ((vClmb(j) < 0) & (clmbStrt > 0))
     kC1 = kC1 + 1;
     mIndx(kC1)= j;
   end    
end

% convert from m/s to km/h
vHn = 3.6*vHn;



hgtMaxHvr = hgt(1);
vHn2 = vHn;
vHn(mIndx)   = [];

vHn2(mIndx2) = [];

vClmb(mIndx) = [];
vDsc(mIndx)  = [];

hgt(mIndx2)   = [];
cTH(mIndx2)   = [];
pReq(mIndx2)  = [];
qReq= pReq/oM;

if vMaxN > 0
   fprintf('vMaxN %5.1f [km/h]    %5.1f [mph]\n',vMaxN*3.6,vMaxN*3.6/1.609);
else
   fprintf('ERROR: vMaxN not > 0. Probably insufficient power\n');
end


fprintf('hover Height  OGE %6.1f [m]   %6.1f [ft]\n',hgtMaxHvr,hgtMaxHvr/f2m);
fprintf('hover Height  IGE %6.1f [m]   %6.1f [ft]\n',hgtHvrIG,hgtHvrIG/f2m);
 
if vMaxH > 0
  fprintf('vMaxH  %5.1f [km/h]   at height  %5.1f [m]\n',vMaxH*3.6,hMaxV);
  fprintf('vMaxH  %5.1f [mph]    at height  %5.1f [ft]\n',vMaxH*3.6/1.609,hMaxV/f2m);
  fprintf('service ceiling  %5.1f [m]    %5.1f [ft]\n',hgtMax,hgtMax/f2m);
end



figure(3)
plot(vHn, vClmb);
xlabel('[km/h]');
ylabel('[m/s]');


outStr = ['Rate of Climb vs. Speed   '  sprintf('H= %5.0f feet',hFlight/f2m)];
title(outStr)

grid

figure(4)
plot(vHn, vDsc);
xlabel('[km/h]');
ylabel('[m/s]');
outStr = ['Rate of Descent vs. Speed   '  sprintf('H= %5.0f feet',hFlight/f2m)];
title(outStr)
grid

figure(5)
title('Flight Envelope')
plot(vHn2, hgt);
grid

figure(6)
title('Thrust Coefficient')
plot(vHn2, cTH);
grid

figure(7)
title('Rotor Torque')
plot(vHn2, qReq);
grid
