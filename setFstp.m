  function j1 = setFstp(jStp,tOffs,sc)
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
% This modul implements varius shapes of control functions, e.g. a step 
% function or saw tooth, for use in the Simulink model of RotaryWingSim
%
% Input
%  jStp  number of function type to be used
%  tOffs time offset at which control input to model starts
%  sc    scale factor for control amplitude
% 
% Output
%   j1   number of function type used  
% 
  global fStp1
  j1=jStp;
  switch jStp
    case {1}
      % uneaven saw tooth
      fStp1(1,1)=0.0; fStp1(1,2)=2.75; fStp1(1,3)=3.0; fStp1(1,4)= 4.5; fStp1(1,5)= 5.0; fStp1(1,6)= 5.5;  fStp1(1,7)=7.0; fStp1(1,8)= 7.25;  fStp1(1,9)= 7.5; fStp1(1,10)= 8.0; fStp1(1,11)= 8.5; fStp1(1,12)= 9.35; fStp1(1,13)= 9.65; fStp1(1,14)= 120;
      fStp1(2,1)=0.0; fStp1(2,2)=0.0; fStp1(2,3)=-5.2; fStp1(2,4)=-4.0; fStp1(2,5)= 5.5; fStp1(2,6)= 5.8;  fStp1(2,7)=0.0; fStp1(2,8)= -7.0;  fStp1(2,9)=-7.8; fStp1(2,10)=-6.5; fStp1(2,11)= 1.5; fStp1(2,12)= 1.25; fStp1(2,13)=  0.0; fStp1(2,14)= 0.0;
    case {2}
      % saw tooth
      fStp1(1,1)=0.0; fStp1(1,2)=3.0;  fStp1(1,3)=3.745; fStp1(1,4)=3.75; fStp1(1,5)=5.248; fStp1(1,6)=5.25;  fStp1(1,7)=6.0; fStp1(1,8)=500;   
      fStp1(2,1)=0.0; fStp1(2,2)=0.0;  fStp1(2,3)=1.0; fStp1(2,4)=1.0; fStp1(2,5)=-1.0; fStp1(2,6)=-1.0;  fStp1(2,7)=0.0; fStp1(2,8)=0.0;
    case {3}
      % [0 3  3.0001    8.0   8.000001       13.0    13.00001                500]
      % [0 0  1.0000   1.0      -1.0         -1.0         0                    0]
      fStp1(1,1)=0.0; fStp1(1,2)=3.0;  fStp1(1,3)=3.001; fStp1(1,4)=8.0; fStp1(1,5)=8.001; fStp1(1,6)=10;  fStp1(1,7)=20.0; fStp1(1,8)=500;   
      fStp1(2,1)=0.0; fStp1(2,2)=0.0;  fStp1(2,3)=1.0; fStp1(2,4)=1.0; fStp1(2,5)=0.0; fStp1(2,6)=0.0;  fStp1(2,7)=0.0; fStp1(2,8)=0.0;
    case {4}
      % saw tooth Montgomerie Parsons Glasgow p. 121
      t1= 0.3; t2=0.5;
      fStp1(1,1)=0.0; fStp1(1,2)=1.15; fStp1(1,3)=1.7; fStp1(1,4)=1.725; fStp1(1,5)= 2.6; fStp1(1,6)= 2.75;  fStp1(1,7)=3.20;  fStp1(1,8)= 4.3 + t2;     fStp1(1,9)=   8 + t2;   fStp1(1,10)= 15;   fStp1(1,11)=500;   
      fStp1(2,1)=0.0; fStp1(2,2)=0.0;  fStp1(2,3)=1.0; fStp1(2,4)=1.0;   fStp1(2,5)=-1.8; fStp1(2,6)=-1.8;   fStp1(2,7)=0.15+t1;   fStp1(2,8)=0.15+ t1;    fStp1(2,9)=-0.1;   fStp1(2,10)=0.0;   fStp1(2,11)=0;   

    case {5}
      % saw tooth Montgomerie Parsons Glasgow p. 121
      t1= 0.3; t2=0.0;
      fStp1(1,1)=0.0; fStp1(1,2)=1.15; fStp1(1,3)=1.7; fStp1(1,4)=1.725; fStp1(1,5)= 2.6; fStp1(1,6)= 2.75;  fStp1(1,7)=3.20;  fStp1(1,8)= 4.0;         fStp1(1,9)= 8;   fStp1(1,10)= 15;   fStp1(1,11)=500;   
      fStp1(2,1)=0.0; fStp1(2,2)=0.0;  fStp1(2,3)=1.0; fStp1(2,4)=1.0;   fStp1(2,5)=-2.4; fStp1(2,6)=-2.4;   fStp1(2,7)=0.15+t1;   fStp1(2,8)=0.15+t1;  fStp1(2,9)= 0;   fStp1(2,10)=0.0;   fStp1(2,11)=0;   
      
  otherwise
      % 5 sec step input
      fStp1(1,1)=0.0; fStp1(1,2)=2.9999;  fStp1(1,3)=3.0; fStp1(1,4)=7.99999; fStp1(1,5)=8.0; fStp1(1,6)=10.0;  fStp1(1,7)=16.0; fStp1(1,8)=500;
      fStp1(2,1)=0.0; fStp1(2,2)=0.0   ;  fStp1(2,3)=1.0; fStp1(2,4)=1.0;     fStp1(2,6)=0.0; fStp1(2,6)=0.0;  fStp1(2,7)=0.0; fStp1(2,8)=0.0;  
  end    

  
  for j=2:length(fStp1)
          fStp1(1,j)=fStp1(1,j)+tOffs;
  end
 
  fStp1(2,:)=sc*fStp1(2,:);
 