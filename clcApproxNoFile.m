function clcApproxNoFile(Xu, Xw, Xq, Zu, Zw, Zq, Mu, Mw, Mq, vBdy)
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
%         --- Version 2.0 ---
%
% This modul calculates approximations for the aircraft modes
% and also the mode values using the full 3x3 stabiliy matrix
% Input:
%  all stability derivatives
%
% Output:
%  to file
%
global gEarth mHeli Iyy 
 
  uB=vBdy(1); vB=vBdy(2); wB=vBdy(3);
  
%   Xu=bMd(1,1);  Xw=bMd(1,2);  Xq=bMd(1,3); 
%   Zu=bMd(2,1);  Zw=bMd(2,2);  Zq=bMd(2,3);
%   Mu=bMd(3,1);  Mw=bMd(3,2);  Mq=bMd(3,3); 


Xq=Xq - wB;  Zq= Zq+ uB;
  
  Zwd=0;  Mwd=0;  gE=gEarth;
  
 
jTyp=2;
if jTyp==1
	fprintf('Xu     %10.6f    %10.2f\n',Xu,Xu*mHeli);
	fprintf('Xw     %10.6f    %10.2f\n',Xw,Xw*mHeli);  
	fprintf('Xq     %10.6f    %10.2f\n',Xq+wB,(Xq+wB)*Iyy);  
                                  
                                  
	fprintf('Zu     %10.6f    %10.2f\n',Zu,Zu*mHeli);
	fprintf('Zw     %10.6f    %10.2f\n',Zw,Zw*mHeli);
	fprintf('Zwd    %10.6f    %10.2f\n',Zwd,Zwd*mHeli);
	fprintf('Zq     %10.6f    %10.2f\n',Zq-uB,(Zq-uB)*Iyy);
	                                
	fprintf('Mu     %10.6f    %10.2f\n',Mu,Mu*mHeli);
	fprintf('Mw     %10.6f    %10.2f\n',Mw,Mw*mHeli);
	fprintf('Mwd    %10.6f    %10.2f\n',Mwd,Mwd*mHeli);
	fprintf('Mq     %10.6f    %10.2f\n',Mq,Mq*Iyy);
else
	fprintf('Xu     %10.2f    %10.8f\n',Xu,Xu/mHeli);
	fprintf('Xw     %10.2f    %10.8f\n',Xw,Xw/mHeli);  
	fprintf('Xq     %10.2f    %10.8f\n',Xq+wB,(Xq+wB)/Iyy);  
                                  
                                  
	fprintf('Zu     %10.2f    %10.8f\n',Zu,Zu/mHeli);
	fprintf('Zw     %10.2f    %10.8f\n',Zw,Zw/mHeli);
	fprintf('Zwd    %10.2f    %10.8f\n',Zwd,Zwd/mHeli);
	fprintf('Zq     %10.2f    %10.8f\n',Zq-uB,(Zq-uB)/Iyy);
	                                
	fprintf('Mu     %10.2f    %10.8f\n',Mu,Mu/mHeli);
	fprintf('Mw     %10.2f    %10.8f\n',Mw,Mw/mHeli);
	fprintf('Mwd    %10.2f    %10.8f\n',Mwd,Mwd/mHeli);
	fprintf('Mq     %10.2f    %10.8f\n',Mq,Mq/Iyy);
    
end
xh=Zu*Mw - Zw*Mu;
fprintf('E = %9.5f   %9.2f\n',xh,xh*mHeli^2);
% from "Padfield, Helicopter Flight Dynamics"
fprintf('\n##### SP and PH Mode ----- approximation ###\n');
fprintf('----- Short Period Mode -----\n');

% Zq below is actually Zq + uB
oMg2_SP  = Mq*Zw - Zq*Mw;
% the positive signs below have been changed 
% to match results from quartic
zT2oM_SP = -Zw-Mq;

zOut=1;

if (oMg2_SP < 0)
  fprintf('negative oMg2_SP\n');
  zOut=0;
else
  oMg_SP=sqrt(oMg2_SP);
  fprintf('%s %12.4f\n','omega SP  ',oMg_SP);
  zt_SP=0.5*zT2oM_SP/oMg_SP;
  fprintf('%s %12.4f\n','zeta  SP  ',zt_SP);

	xh= 0.25*zT2oM_SP^2-oMg2_SP;
	
	if (xh < 0)
	  % fprintf('negative root\n');
	  xh=-xh;
	  reSP = -0.5*zT2oM_SP; imSP=sqrt(xh);
	  fprintf('re %8.4f +%8.4fj\n',reSP,imSP);
	  fprintf('re %8.4f -%8.4fj\n',reSP,imSP);
	else
	  reSP = -0.5*zT2oM_SP; imSP=sqrt(xh);
	  fprintf('re %8.4f + %8.4f\n',reSP,imSP);
	  fprintf('re %8.4f - %8.4f\n',reSP,imSP);
	end
end  



fprintf('\n----- Phugoid Mode -----\n');
% Zq below is actually Zq + uB
oMg2_PH  = -gEarth/uB*(Zu - Zw*(Zu*Mq - Mu*Zq)/oMg2_SP);
% Xq below is actually Xq - wB
zT2oM_PH = -Xu + ((Xw - gEarth/uB)*(Zu*Mq - Mu*Zq) + Xq*(Zw*Mu - Mw*Zu))/oMg2_SP;


if (oMg2_PH < 0)
  fprintf('negative oMg2_PH\n');
  zOut=0;  
else
  oMg_PH=sqrt(oMg2_PH);
  fprintf('%s %12.4f\n','omega PH  ',oMg_PH);
  zt_PH=0.5*zT2oM_PH/oMg_PH;

  fprintf('%s %12.4f\n','zeta  PH  ',zt_PH);

	xh= 0.25*zT2oM_PH^2-oMg2_PH;
	
	if (xh < 0)
	  % fprintf('negative root\n');
	  xh=-xh;
	  rePH = -0.5*zT2oM_PH; imPH=sqrt(xh);
	  fprintf('re %8.4f +%8.4fj\n',rePH,imPH);
	  fprintf('re %8.4f -%8.4fj\n',rePH,imPH);
	else
	  rePH = -0.5*zT2oM_PH; imPH=sqrt(xh);
	  fprintf('re %8.4f + %8.4f\n',rePH,imPH);
	  fprintf('re %8.4f - %8.4f\n',rePH,imPH);
	end
end  

% syms zS;
if (zOut > 0.5)
z1=rePH + imPH*i;
z2=rePH - imPH*i;
z3=reSP + imSP*i;
z4=reSP - imSP*i;

% f=(zS-z1)*(zS-z2)*(zS-z3)*(zS-z4);
f2(1)= 1;
f2(2)= (-z1-z2-z3-z4);
f2(3)= (z1*z2-(-z1-z2)*z3-(-z1-z2-z3)*z4);
f2(4)= (-z1*z2*z3-(z1*z2-(-z1-z2)*z3)*z4);
f2(5)= z1*z2*z3*z4;

% fc=collect(f,'zS');
% f2=sym2poly(fc);

fprintf('\n');
fprintf('A    %12.4f\n',f2(1));
fprintf('B    %12.4f\n',f2(2));
fprintf('C    %12.4f\n',f2(3));
fprintf('D    %12.4f\n',f2(4));
fprintf('E    %12.4f\n',f2(5));
end

fprintf('\n##### SP and PH Mode from quartic #####\n');
   
cthta =  0.0;
Mthta = -180*cthta;
% 
%        u   w     q    thta
det1 = [ Xu  Xw    Xq   -gE;  % u
         Zu  Zw    Zq    0;  % w
         Mu  Mw    Mq    0;  % q  
         0    0    1     0]; % thta

%      uDot wDot qDot thtaDot  
det2 = [ 1   0   0   0;  % zDot
         0   1   0   0;
         0   0   1   0;    
         0   0   0   1]; % qDot


Xwd=0; Zwd=0; Mwd=0;
[a4 a3 a2 a1 aN]=coeffnc(Xu,Xw,Xq,Xwd,Zu,Zw,Zq,Zwd,Mu,Mw,Mwd,Mq,vBy(1));

%ak(5)=a4;  ak(4)=a3;  ak(3)=a2;  ak(2)=a1;  ak(1)=aN;
ak(1)=a4;  ak(2)=a3;  ak(3)=a2;  ak(4)=a1;  ak(5)=aN;
r=roots(ak);

d2(1)=a4;  d2(2)=a3;  d2(3)=a2;  d2(4)=a1;  d2(5)=aN;


fprintf('A %12.3f  %12.3f\n',d2(1),ak(1));
fprintf('B %12.3f  %12.3f\n',d2(2),ak(2));
fprintf('C %12.3f  %12.3f\n',d2(3),ak(3));
fprintf('D %12.3f  %12.3f\n',d2(4),ak(4));
fprintf('E %12.3f  %12.3f\n\n',d2(5),ak(5));


rthD = d2(2)*d2(3)*d2(4)-d2(1)*d2(4)^2 - d2(2)^2*d2(5);
if abs(rthD) > 5.0e5
  fprintf('%s %12.3f %s\n\n','Routh D',rthD/1.0e+6,'x10^6');
else
  fprintf('%s %12.3f\n\n','Routh D',rthD);
end
fprintf('%s%12.3f \n\n','M thta  ',Mthta);

z=roots(d2);
z1=z(1) + j*0; z2=z(2); z3=z(3); z4=z(4);
fprintf('z1 %9.5f  %9.5fj\n',real(z1),imag(z1));
fprintf('z2 %9.5f  %9.5fj\n',real(z2),imag(z2));
fprintf('z3 %9.5f  %9.5fj\n',real(z3),imag(z3));
fprintf('z4 %9.5f  %9.5fj\n',real(z4),imag(z4));

% fprintf('z1 %9.5f  %9.5fj\n',real(z1),imag(z1));
% fprintf('z2 %9.5f  %9.5fj\n',real(z2),imag(z2));
% fprintf('z3 %9.5f  %9.5fj\n',real(z3),imag(z3));
% fprintf('z4 %9.5f  %9.5fj\n\n',real(z4),imag(z4));
% 
if abs(rthD) > 5.0e5
  fprintf('%s %12.3f %s\n\n','Routh D',rthD/1.0e+6,'x10^6');
else
  fprintf('%s %12.3f\n\n','Routh D',rthD);
end

f2(1)= 1;
f2(2)= (-z1-z2-z3-z4);
f2(3)= (z1*z2-(-z1-z2)*z3-(-z1-z2-z3)*z4);
f2(4)= (-z1*z2*z3-(z1*z2-(-z1-z2)*z3)*z4);
f2(5)= z1*z2*z3*z4;


% fprintf('Backwards\n');
% fprintf('A %12.4f\n',f2(1));
% fprintf('B %12.4f\n',f2(2));
% fprintf('C %12.4f\n',f2(3));
% fprintf('D %12.4f\n',f2(4));
% fprintf('E %12.4f\n',f2(5));
% fprintf('\n');
% 
% 

AL = -inv(det2)*det1;
LongitudinalSystem = ss(AL,[0;0;0;0],eye(4),[0;0;0;0]);
[VL,DL] = eig(AL);
ShortPeriodIC = VL(:,1)+VL(:,2);
PhugoidIC = VL(:,3)+VL(:,4);

fprintf('\n       Short Period                Phugoid\n');
for (j=1:4)
  fprintf('%d %9.5f %8.5fj      %9.5f %8.5fj\n',j,real(ShortPeriodIC(j)),imag(ShortPeriodIC(j)),real(PhugoidIC(j)),imag(PhugoidIC(j)));
end
fprintf('\n');

figure
xh=abs(imag(z1)) + abs(imag(z2));
if abs(xh) < 1.0e-10 
  return
end


[x,t,y] =  initial(LongitudinalSystem,ShortPeriodIC);
plot(x,y)

initial(LongitudinalSystem, PhugoidIC)
title('Pure Phugoid Response')

jImag=-1;
jReal = [-1 -1];
k=0;

% find out where in DL the real solutions are
for (j=1:4)
  if (abs(imag(DL(j,j))) > 1.0e-8)
     jImag=j;
  else
     k=k+1;
     jReal(k)=j;
 end
end

% if there are two complex solutions 
% the larger real value is that of
% the ----- Phugoid Mode -----
if (k==0)
   if (real(DL(1,1)) > real(DL(3,3)))
     jReal(1)=1;
     jImag = 3;
   else
     jReal(1)=3;
     jImag =1;
   end  
   k=1;
end   

DL_PH=DL(jImag,jImag);

fprintf('\n%s\n','----- Phugoid Mode -----');
omegan_PH = abs(DL_PH);
fprintf('%s %12.4f\n','omega PH  ',omegan_PH);
zeta_PH = -real(DL_PH)/omegan_PH;
fprintf('%s %12.4f\n','zeta  PH  ',zeta_PH);


xNN=omegan_PH*sqrt(1-zeta_PH^2);
if (abs(xNN) > 1.0e-8) 
  period_PH = 2*pi/xNN;
  fprintf('%s %12.4f\n','period PH ',period_PH);
else
  fprintf('%s\n','period PH   --------');
end


t_half_PH = abs(log(0.5)/real(DL_PH));
fprintf('%s %12.4f\n','tHalf  PH ',t_half_PH);
N_half_PH = abs((log(0.5)/(2*pi))*(imag(DL_PH)/real(DL_PH)));
fprintf('%s %12.4f\n','N half PH ',N_half_PH);

fprintf ('\n%s\n','----- Short Period Mode -----');
for (i=1:k)
DL_SP=DL(jReal(i),jReal(i));

omegan_SP = abs(DL_SP);
fprintf('%s %12.4f\n','omega SP  ',omegan_SP);
zeta_SP = -real(DL_SP)/omegan_SP;
fprintf('%s %12.4f\n','zeta  SP  ',zeta_SP);
xNN=omegan_SP*sqrt(1-zeta_SP^2);
if (abs(xNN) > 1.0e-6) 
  period_SP = 2*pi/xNN;
  fprintf('%s %12.4f\n','period SP ',period_SP);
else
  fprintf('%s\n','period SP   --------');
end

t_half_SP = abs(log(0.5)/real(DL_SP));
fprintf('%s %12.4f\n','tHalf  SP ',t_half_SP);
N_half_SP = abs((log(0.5)/(2*pi))*(imag(DL_SP)/real(DL_SP)));
fprintf('%s %12.4f\n\n','N half SP ',N_half_SP);

end
