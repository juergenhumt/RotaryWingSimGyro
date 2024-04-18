function readFStp(inFileName,tConst,tEnd,aMax)
global fStp1
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
% this module reads control input vs time data to enable control input sequences 
% as measured in flight. Data are stored in global variable fStp1
%
% Input
%  inFileName  input file name
%  tConst      time offset value for which controls remain constant
%  tEnd        time value of end of control input
%  aMax        amplitude scaling factor 
% 
%
fAux=csvread(inFileName);

iEnd=length(fAux);


retMin=min(fAux,[],1);
retMax=max(fAux,[],1);



iStrt=0;
tEnd=tEnd-tConst;


if tConst > 0
  fStp1(1,1)=0;
  fStp1(1,2)=0;
  iStrt=1;
end

tScale=tEnd/(retMax(1)-retMin(1));
aScale= max(abs(retMax(2)),abs(retMin(2)));
aScale=aMax/(aScale-fAux(1,2));

for i=1:iEnd
   i1=i+iStrt;
   
   fStp1(1,i1)=tScale*(fAux(i,1)-retMin(1)) + tConst;
   fStp1(2,i1)=aScale*(fAux(i,2)-fAux(1,2));
end

figure(5)
plot(fStp1(1,:),fStp1(2,:))
