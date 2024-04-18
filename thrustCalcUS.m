function tProp =thrustCalcUS(vGyroInp, thrtlStng)
%
% Copyright 2011 Juergen Humt
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
% Input:
%   vFlight   : flight speed [m/s] converted to [miles/h]
%   thrtlStng : throttle setting, mulitplying factor for engine power
%               0 <= thrtlStn <= 1.0 
%
% Global:
%   engineRPM : engine speed [rpm]
%   dProp     : propeller diameter [ft]
%
% Output:
%   tProp : propeller thrust in Newton!! [N] 
%

global f2m slg2kg lbMass lbForce rhoAir gEarth pCTP propDia engineRPM engineHP

cPr = 3.325e10*engineHP/(engineRPM^3*propDia^5);
cTP= polyval(pCTP,cPr);

tS = cTP*thrtlStng*engineHP*33000*0.9/(engineRPM*propDia);

% speed in the propeller thrust equation is in miles per hour
vFlight=vGyroInp*3.6/1.605;

jCP=1.237*vFlight*sqrt(engineRPM*propDia^3/(engineHP*1.0e7));

aTTS=-0.4/3;  bTTS = 1.0;
tTS = aTTS*jCP+bTTS;


tProp = tTS*tS*lbForce;
