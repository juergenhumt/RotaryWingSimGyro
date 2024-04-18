function clcApprox(Xu, Xw, Xq, Zu, Zw, Zq, Mu, Mw, Mq, vBdy, mu, alfaD, gam, cTs, outFile)
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
global gEarth mHeli Iyy octRun f2m

cAux=1.35;

uB=vBdy(1); vB=vBdy(2); wB=vBdy(3);
  

Xq=Xq - wB;  Zq= Zq+ uB;
  
  Zwd=0;  Mwd=0;  gE=gEarth;
    

	fprintf(outFile,'Xu     %10.6f    %10.2f    %10.2f\n',Xu,Xu*mHeli,Xu*mHeli*f2m);
	fprintf(outFile,'Xw     %10.6f    %10.2f    %10.2f\n',Xw,Xw*mHeli,Xw*mHeli*f2m);  
	fprintf(outFile,'Xq     %10.6f    %10.2f    %10.2f\n',Xq+wB,(Xq+wB)*Iyy,(Xq+wB)*Iyy*cAux);  
                                  
                                  
	fprintf(outFile,'Zu     %10.6f    %10.2f    %10.2f\n',Zu,Zu*mHeli,Zu*mHeli*f2m);
	fprintf(outFile,'Zw     %10.6f    %10.2f    %10.2f\n',Zw,Zw*mHeli,Zw*mHeli*f2m);
	fprintf(outFile,'Zwd    %10.6f    %10.2f    %10.2f\n',Zwd,Zwd*mHeli,Zwd*mHeli*f2m);
	fprintf(outFile,'Zq     %10.6f    %10.2f    %10.2f\n',Zq-uB,(Zq-uB)*Iyy,(Zq-uB)*Iyy*cAux);
	                                
	fprintf(outFile,'Mu     %10.6f    %10.2f    %10.2f\n',Mu,Mu*mHeli,Mu*mHeli*f2m);
	fprintf(outFile,'Mw     %10.6f    %10.2f    %10.2f\n',Mw,Mw*mHeli,Mw*mHeli*f2m);
	fprintf(outFile,'Mwd    %10.6f    %10.2f    %10.2f\n',Mwd,Mwd*mHeli,Mwd*mHeli*f2m);
	fprintf(outFile,'Mq     %10.6f    %10.2f    %10.2f\n',Mq,Mq*Iyy,Mq*Iyy*cAux);
	
xh=Zu*Mw - Zw*Mu;
fprintf('E = %9.5f   %9.2f\n',xh,xh*mHeli^2);
% from "Padfield, Helicopter Flight Dynamics"
fprintf(outFile,'\n##### SP and PH Mode ----- approximation ###\n');
fprintf(outFile,'----- Short Period Mode -----\n');


% Zq below is actually Zq + uB
oMg2_SP  = Zw*Mq - Zq*Mw;  % Padfield 4.144, p 243
% the positive signs below have been changed 
% to match results from quartic
zT2oM_SP = -Zw-Mq;    % Padfield 4.144, p 243

zOut=1;

if (oMg2_SP < 0)
  fprintf(outFile,'negative oMg2_SP\n');
  zOut=0;
else
  oMg_SP=sqrt(oMg2_SP);
  fprintf(outFile,'%s %12.4f\n','omega SP  ',oMg_SP);
  zt_SP=0.5*zT2oM_SP/oMg_SP;
  fprintf(outFile,'%s %12.4f\n','zeta  SP  ',zt_SP);

	xh= 0.25*zT2oM_SP^2-oMg2_SP;
	
	if (xh < 0)
	  % fprintf(outFile,'negative root\n');
	  xh=-xh;
  end
	  reSP = -0.5*zT2oM_SP; imSP=sqrt(xh);
	  fprintf(outFile,'re %8.4f +%8.4fj\n',reSP,imSP);
	  fprintf(outFile,'re %8.4f -%8.4fj\n',reSP,imSP);
end  



fprintf(outFile,'\n----- Phugoid Mode -----\n');
% Zq below is actually Zq + uB
oMg2_PH  = -gEarth/uB*(Zu - Zw*(Zu*Mq - Mu*Zq)/oMg2_SP);
% Xq below is actually Xq - wB
zT2oM_PH = -Xu + ((Xw - gEarth/uB)*(Zu*Mq - Mu*Zq) + Xq*(Zw*Mu - Mw*Zu))/oMg2_SP;


if (oMg2_PH < 0)
  fprintf(outFile,'negative oMg2_PH\n');
  zOut=0;  
else
  oMg_PH=sqrt(oMg2_PH);
  fprintf(outFile,'%s %12.4f\n','omega PH  ',oMg_PH);
  zt_PH=0.5*zT2oM_PH/oMg_PH;

  fprintf(outFile,'%s %12.4f\n','zeta  PH  ',zt_PH);

	xh= 0.25*zT2oM_PH^2-oMg2_PH;
	
	if (xh < 0)
	  % fprintf(outFile,'negative root\n');
	  xh=-xh;
	  rePH = -0.5*zT2oM_PH; imPH=sqrt(xh);
	  fprintf(outFile,'re %8.4f +%8.4fj\n',rePH,imPH);
	  fprintf(outFile,'re %8.4f -%8.4fj\n',rePH,imPH);
	else
	  rePH = -0.5*zT2oM_PH; imPH=sqrt(xh);
	  fprintf(outFile,'re %8.4f + %8.4f\n',rePH,imPH);
	  fprintf(outFile,'re %8.4f - %8.4f\n',rePH,imPH);
	end
end  

fprintf(outFile,'\n##### SP and PH Mode from quartic 3104 #####\n');
   
Xwd=0; Zwd=0; Mwd=0; Mthta=0;
[A B C D E]=coeffnc3104 (Xu,Xw,Xq,Xwd,Zu,Zw,Zq,Zwd,Mu,Mw,Mq,Mwd,Mthta,mu,alfaD,gam,cTs);

d2(1)=A;  d2(2)=B;  d2(3)=C;  d2(4)=D;  d2(5)= E;

fprintf(outFile,'quartic coefficients\n');
fprintf(outFile,'A %12.3f\n',d2(1));
fprintf(outFile,'B %12.3f\n',d2(2));
fprintf(outFile,'C %12.3f\n',d2(3));
fprintf(outFile,'D %12.3f\n',d2(4));
fprintf(outFile,'E %12.3f\n\n',d2(5));


rthD = d2(2)*d2(3)*d2(4)-d2(1)*d2(4)^2 - d2(2)^2*d2(5);
if abs(rthD) > 5.0e5
  fprintf(outFile,'%s %12.3f %s\n\n','Routh D',rthD/1.0e+6,'x10^6');
else
  fprintf(outFile,'%s %12.3f\n\n','Routh D',rthD);
end
fprintf(outFile,'%s%12.3f \n\n','M thta  ',Mthta);

z=roots(d2);
z1=z(1); z2=z(2); z3=z(3); z4=z(4);
fprintf(outFile,'z1 %9.5f  %9.5fj\n',real(z1),imag(z1));
fprintf(outFile,'z2 %9.5f  %9.5fj\n',real(z2),imag(z2));
fprintf(outFile,'z3 %9.5f  %9.5fj\n',real(z3),imag(z3));
fprintf(outFile,'z4 %9.5f  %9.5fj\n',real(z4),imag(z4));

fprintf('z1 %9.5f  %9.5fj\n',real(z1),imag(z1));
fprintf('z2 %9.5f  %9.5fj\n',real(z2),imag(z2));
fprintf('z3 %9.5f  %9.5fj\n',real(z3),imag(z3));
fprintf('z4 %9.5f  %9.5fj\n\n',real(z4),imag(z4));

if abs(rthD) > 5.0e5
  fprintf('%s %12.3f %s\n\n','Routh D',rthD/1.0e+6,'x10^6');
else
  fprintf('%s %12.3f\n\n','Routh D',rthD);
end


AL =[Xu                     Xw*sin(gam)              Xq               -gEarth*cos(gam);
    [Zu                     Zw*cos(gam)              Zq     -gEarth*sin(gam)]/(1-Zwd);
    [Mu+Zu*Mwd/(1-Zwd) Mw+Zw*Mwd/(1-Zwd) Mq+Zq*Mwd/(1-Zwd) ...
         -gEarth*sin(gam)*Mwd/(1-Zwd)];
       [ 0 0 1 0]];

damp(AL)

LongitudinalSystem = ss(AL,[0;0;0;0],eye(4),[0;0;0;0]);
[VL,DL] = eig(AL);
ShortPeriodIC = VL(:,1)+VL(:,2);
PhugoidIC = VL(:,3)+VL(:,4);

fprintf(outFile,'\n       Short Period                Phugoid\n');

for j=1:4
fprintf('%8.5f   ',real(DL(j,j)));
fprintf(outFile,'%8.5f   ',real(DL(j,j)));

if abs(imag(DL(j,j))) > 1.0e-12
  fprintf('%8.5fj', imag(DL(j,j)));
  fprintf(outFile,'%8.5fj', imag(DL(j,j)));  
end
fprintf('\n');
fprintf(outFile,'\n');
end
fprintf('\n');
fprintf(outFile,'\n');


for j=1:4
fprintf('%8.5f   ',real(z(j)));
fprintf(outFile,'%8.5f   ',real(z(j)));

if abs(imag(z(j))) > 1.0e-12
  fprintf('%8.5fj', imag(z(j)));
  fprintf(outFile,'%8.5fj', imag(z(j)));  
end
fprintf('\n');
fprintf(outFile,'\n');
end

ssp1=(-B+sqrt(B*B-4*A*C))/2/A;
ssp2=(-B-sqrt(B*B-4*A*C))/2/A;


h1=figure; % (jFig);
t=[0:0.05:10];
d=real(exp(ssp1*t));
plot(t,d);
title('Short Period');
xlabel('time (sec)');
ylabel('unit disturbance');
text(2,0.6,'s = ');
output=num2str(ssp1);
text(3,0.6,output);

A1=1.0;
B1=D/C - E*B/C/C;
C1=E/C;

ph1=(-B1+sqrt(B1*B1-4*A1*C1))/2/A1;
ph2=(-B1-sqrt(B1*B1-4*A1*C1))/2/A1;

% jFig=jFig+1;
h2=figure; % (jFig);
t=[0:1:1000];
d=real(exp(ph1*t));
plot(t,d);
title('Phugoid');
xlabel('time (sec)');
ylabel('unit disturbance');
text(200,0.6,'s = ');
output=num2str(ph1);
text(300,0.6,output);

% jFig=jFig+1;
h3=figure; % (jFig);
t=[0:0.01:100];
d=real(exp(ph1*t))+real(exp(ssp1*t));
plot(t,d);
title('Combined Motion');
xlabel('time (sec)');
ylabel('unit disturbance');
% ==============================================


for (j=1:4)
  fprintf(outFile,'%d %9.5f %8.5fj      %9.5f %8.5fj\n',j,real(ShortPeriodIC(j)),imag(ShortPeriodIC(j)),real(PhugoidIC(j)),imag(PhugoidIC(j)));
end
fprintf(outFile,'\n');

figure
% xh=abs(imag(z1)) + abs(imag(z2));
% if abs(xh) < 1.0e-10 
%   return
% end

% if octRun > 0.5
%   [x,t,y] =  initial(LongitudinalSystem,ShortPeriodIC);
%   % initial(LongitudinalSystem,ShortPeriodIC)
%   plot(x,y)
% else
%     initial(LongitudinalSystem,ShortPeriodIC)
% end
% 
% initial(LongitudinalSystem, PhugoidIC)
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

fprintf(outFile,'\n%s\n','----- Phugoid Mode -----');
omegan_PH = abs(DL_PH);
fprintf(outFile,'%s %12.4f\n','omega PH  ',omegan_PH);
zeta_PH = -real(DL_PH)/omegan_PH;
fprintf(outFile,'%s %12.4f\n','zeta  PH  ',zeta_PH);
fhi_PH = rad2deg(atan(sqrt(1 - zeta_PH^2)/zeta_PH));
fprintf(outFile,'%s %12.4f\n','fhi  PH  ',fhi_PH);



xNN=omegan_PH*sqrt(1-zeta_PH^2);
if (abs(xNN) > 1.0e-8) 
  period_PH = 2*pi/xNN;
  fprintf(outFile,'%s %12.4f\n','period PH ',period_PH);
else
  fprintf(outFile,'%s\n','period PH   --------');
end



t_half_PH = abs(log(0.5)/real(DL_PH));
fprintf(outFile,'%s %12.4f\n','tHalf  PH ',t_half_PH);
N_half_PH = abs((log(0.5)/(2*pi))*(imag(DL_PH)/real(DL_PH)));
fprintf(outFile,'%s %12.4f\n','N half PH ',N_half_PH);
