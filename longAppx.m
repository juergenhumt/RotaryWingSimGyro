function longAppx(Xu, Xw, Xq, Zu, Zw, Zq, Mu, Mw, Mq, vBdy, jFig)
global gEarth

Vx = vBdy(1);
Mwdot=0;

A=1.0;
B=-Mq-Mwdot*Vx-Zw-Xu;
C=Zw*Mq-Vx*Mw-Xw*Zu+Xu*Mq+Xu*Mwdot*Vx+Xu*Zw;
D=-Xu*Zw*Mq+Xu*Mw*Vx+Zu*Xw*Mq+Zu*gEarth*Mwdot-Mu*Xw*Vx+Mu*gEarth;
E=gEarth*Zu*Mw-gEarth*Mu*Zw;

Xde = 0;
Zde = 92.0;
Mde = -7.10;

At = Mde + Zde*Mwdot;
Bt = Xde*(Zw*Mwdot + Mu) + Zde*(Mw - Xu*Mwdot) - Mde*(Xu + Zw);
Ct = Xde*(Zu*Mw - Zw*Mu) + Zde*(Mu*Xw - Mw*Xu) + Mde*(Zw*Xu - Xw*Zu);


fprintf('\n');
fprintf('A    %12.4f\n',A);
fprintf('B    %12.4f\n',B);
fprintf('C    %12.4f\n',C);
fprintf('D    %12.4f\n',D);
fprintf('E    %12.4f\n',E);

Xwd=0; Zwd=0; Mwd=0; Mthta=0;
[a4 a3 a2 a1 aN]=coeffnc(Xu,Xw,Xq,Xwd,Zu,Zw,Zq,Zwd,Mu,Mw,Mwd,Mq,Mthta);
% [a4 a3 a2 a1 aN]=coeffnc(Xu,Xw,Xq,Xwd,Zu,Zw,Zq,Zwd,Mu,Mw,Mwd,Mq,Mthta,vBdy(1));

d2(1)=a4;  d2(2)=a3;  d2(3)=a2;  d2(4)=a1;  d2(5)=aN;
r=roots(d2);

fprintf('xxxx----xxxx\n');
fprintf('A    %12.4f\n',d2(1));
fprintf('B    %12.4f\n',d2(2));
fprintf('C    %12.4f\n',d2(3));
fprintf('D    %12.4f\n',d2(4));
fprintf('E    %12.4f\n\n',d2(5));

rthD = d2(2)*d2(3)*d2(4)-d2(1)*d2(4)^2 - d2(2)^2*d2(5);
if abs(rthD) > 5.0e5
  fprintf('%s %12.3f %s\n\n','Routh D',rthD/1.0e+6,'x10^6');
else
  fprintf('%s %12.3f\n\n','Routh D',rthD);
end

P=[A B C D E];
s=roots(P);

for j=1:4
fprintf('%8.5f   ',real(s(j)));
if abs(imag(s(j))) > 1.0e-12
  fprintf('%8.5fj', imag(s(j)));
end
fprintf('\n');
end

ssp1=(-B+sqrt(B*B-4*A*C))/2/A;
ssp2=(-B-sqrt(B*B-4*A*C))/2/A;


h1=figure(jFig);
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

jFig=jFig+1;
h2=figure(jFig);
t=[0:1:1000];
d=real(exp(ph1*t));
plot(t,d);
title('Phugoid');
xlabel('time (sec)');
ylabel('unit disturbance');
text(200,0.6,'s = ');
output=num2str(ph1);
text(300,0.6,output);

jFig=jFig+1;
h3=figure(jFig);
t=[0:0.01:100];
d=real(exp(ph1*t))+real(exp(ssp1*t));
plot(t,d);
title('Combined Motion');
xlabel('time (sec)');
ylabel('unit disturbance');


