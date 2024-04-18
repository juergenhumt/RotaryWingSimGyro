function [a4,a3,a2,a1,aN]= coeffnc(Xu,Xw,Xq,Xwd,Zu,Zw,Zq,Zwd,Mu,Mw,Mwd,Mq,Mthta,Un)
%                                 (Xu,Xw,Xq,Xwd,Zu,Zw,Zq,Zwd,Mu,Mw,Mwd,Mq,Mthta,uB);
global Iyy m gEarth


% a4 = Iyy^2
a4= Iyy^2;


% a3=-Mq-Zw*Iyy^2-Xu*Iyy^2;
a3 = (-Mq-Xu*Iyy^2-Zw*Iyy^2);


% a2=-Mthta+Zw*Mq-Zu*Xw*Iyy^2+Mw*Un+Xu*Mq+Xu*Zw*Iyy^2;
a2 = (Mw*Un-Zu*Xw*Iyy^2+Xu*Zw*Iyy^2-Mthta+Xu*Mq+Zw*Mq);

% a1=-Xu*Zw*Mq-Xu*Zw*Mthta-Xu*Mw*Un+Zw*Mthta+Xu*Mthta+Zu*Xw*Mq;
a1 = (Zw*Mthta+Zu*Xw*Mq+Xu*Mthta-Xu*Zw*Mq-Xu*Mw*Un);

% aN=+Zu*Xw*Mthta+Zu*Mw*gEarth;
aN= Zu*Mw*gEarth+Zu*Xw*Mthta-Xu*Zw*Mthta;



                                                                                                                                                                                                  