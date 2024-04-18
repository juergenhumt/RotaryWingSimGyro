  clc
  close all
  
  
%   gamRot=  10.5;
%   cMIN  =   2.4;
%   thtaN =   2.5; 
%   Gmod  =  800; % 2305.0
%   cGA   =   0.022;


mu=0.2;

k= 3;


if k == 0
% data from report naca 600
  aLift  =    5.73;      
  gamRot =   12.5;
  cMIN   =    2.4;
  thtaN  =    5.5; 
  cGA    =    0.0116;
  CmPrfl =   -0.056;
  Gmod   = 2238.0;
  dDrgN= 0.0085; dDrg1= -0.0416; dDrg2=0.03329;
end 



if k == 1
  aLift  =    5.73;  
  gamRot =   12.5;
  cMIN   =    2.4;
  thtaN  =    6.1;  % KD-1 thtan=5.5°
  Gmod   = 1515.0;
  cGA    =    0.0455;
  CmPrfl =   -0.09;
end


if k == 2
% numbers modified to give a better agreement with 
% measured data of aN, a1, b1
  aLift  =    5.73;      
  gamRot =   12.5;
  cMIN   =    2.4;
  thtaN  =    5.5; 
  cGA    =    0.011;
  CmPrfl =   -0.062;
  Gmod   = 2005.0;
  dDrgN= 0.0165; dDrg1= -0.0416; dDrg2=0.03329;
end 



if k == 3
% numbers modified to give a better agreement with 
% measured data of aN, a1, b1
  aLift  =    5.73;      
  gamRot =   12.5;
  cMIN   =    2.4;
  thtaN  =    5.5; 
  cGA    =    0.011;
  CmPrfl =   -0.07;
  Gmod   = 2305.0;
  dDrgN= 0.0145; dDrg1= -0.0416; dDrg2=0.03329;
end 


% ===============
%   aLift  =    5.73;      
%   gamRot =   12.5;
%   cMIN   =    2.4;
%   thtaN  =    5.5; 
%   cGA    =    0.0116;
%   CmPrfl =   -0.06;
%   Gmod   = 2205.0;
%   dDrgN= 0.0085; dDrg1= -0.0416; dDrg2=0.03329;
% ===============

  c1=1.605/3.6/rRot*30/pi;
  thtaN  = deg2rad(thtaN); 
  
  vGry=[30:1.0:95];
  for j=1:length(vGry)
    rrpmAux(j)= crvVal(vGry(j),rrpmCff);
    muAux(j)  = vGry(j)/rrpmAux(j)*c1;
  end
  
  rrpmAxCff  = polyfit(muAux,rrpmAux,3);
  
  muList=0;
  
  muList=[0.125:0.0125:0.325];
   
   wC=0;

   zB=0;
   pC=0;
   
  for j=1:length(muList)
    mu=muList(j);
   
    aN600(j)= crvVal(mu, aNCff);
    a1600(j)= crvVal(mu, a1Cff);
    b1600(j)= crvVal(mu, b1Cff);
    la600(j)= crvVal(mu, laCff);

   
    % rrpm600(j)= crvVal(mu,rrpmCff);
    rrpm600(j)= polyval(rrpmAxCff, mu);
    vMps    = mu*pi*rrpm600(j)/30*rRot;
    vMph(j) = vMps*3.6/1.605;

    cT600(j)= crvVal(vMph(j),cTCff);
    oM  = pi*rrpm600(j)/30;   
    oMR = oM*rRot;
  
    talfaNf(j) = la600(j)/mu + cT600(j)/(2*mu*sqrt(mu^2 + la600(j)^2)); 
    alfaNf(j) = atan(talfaNf(j));
   

   
    wC = -vMps*sin(alfaNf(j));
  
    [cTs2(j), cT2(j), lamNf(j), lamI(j), kG, jItr]= clcViNf(mu, wC, oM, thtaN, pC, zB, cT600(j));

    TMRx= cT600(j)*(rhoAir*aDiskMR*(oMR)^2);
  
    % [eNTwst(j), e1Twst(j), e2Twst(j), n1Twst(j), n2Twst(j)] = clcCoeffTwst(mu, la600(j), thtaN, oMR, TMRx);
    [eNTwst(j), e1Twst(j), e2Twst(j), n1Twst(j), n2Twst(j)] = clcCoeffTwst(mu, lamNf(j), thtaN, oMR, TMRx);

%     cLa=1.0; ctN= 0.70;  cOM= 1.0; cTM= 1.0;
%     [aN(j), a1r(j), a2r(j), b1r(j), b2r(j), dQp(j)]= clcRotAnglTwst(mu, la600(j)*cLa, thtaN*ctN, oMR*cOM, TMRx*cTM);
%     cLa= 1.0;  ctN=1.0;  cOM= 1.0; cTM= 1.0;
%     [aN(j), a1r(j), a2r(j), b1r(j), b2r(j), dQm(j)]= clcRotAnglTwst(mu, la600(j)*cLa, thtaN*ctN, oMR*cOM, TMRx*cTM);
%     
    [aN(j), a1r(j), a2r(j), b1r(j), b2r(j), dQn(j)]= clcRotAnglTwstKx(mu, la600(j), thtaN, oMR, TMRx);
    cT716(j) = clcCT(mu,la600(j),thtaN,oM,cT600(j));
  

  end
  
  cT716= 0.5*aLift*sigma*cT716;
  
  plot(muList,rrpm600);
  grid
  
  figure(2)
  plot(muList,cT600)
  hold on
  plot(muList,cT2,'g')
  plot(muList,cT716,'y')
  grid
  legend('cT600','cTcalc','cT716')
  
  figure(3)
  plot(muList,la600)
  hold on
  grid
  plot(muList,lamNf,'g')
  legend('lamNf 600','lamNf calc')

%   figure(4)  
%   plot(muList,a1600,'bo-')
  
  
  figure(5)
  plot(muList,rad2deg(eNTwst));
  grid
  hold on
  plot(muList,rad2deg(e1Twst),'go-')
  plot(muList,rad2deg(e2Twst),'gx-')
  plot(muList,rad2deg(n1Twst),'yo-')
  plot(muList,rad2deg(n2Twst),'yx-')

  legend('eN','e1','e2','n1','n2')

  figure(6)
  plot(muList,aN600);
  hold on
  plot(muList,rad2deg(aN),'g')
  plot(muList,a1600,'bo-')
  plot(muList,rad2deg(a1r),'go-')

  plot(muList,b1600,'bx-')
  plot(muList,rad2deg(b1r),'gx-')
  
  figure(7)
  plot(muList,dQn)
  hold on
  grid
%
%   plot(muList,dQm,'g')
%   plot(muList,dQp,'y')

%
% http://www.wwiivehicles.com/unitedkingdom/aircraft/trainer/anson.asp
% http://failheap-challenge.com/showthread.php?892-russian-planes-appreciation-thread/page33
% http://www.youtube.com/watch?v=xyxj8soYqwQ
% http://www.youtube.com/watch?v=4JLUaIEQT0k
% 
% As a further means of increasing simultaneously the lateral, longitudinal and
% directional stability and at the same time economizing in weight and in paraiste
% drag the horizontal tail may, in accordance with this invention, be provided with
% oblique upturned tips
% http://www.youtube.com/watch?v=fjRYdQUTitw&feature=youtu.be
% 
% Miller, Some Aspects Of the Helicopter Stability and Control Problem
%
% K-Max specs:
% Normal rotor rpm is 270, giving maximum blade tip speed of 200m/min; 
% translational lift is attained at 22km/h; rpm reduced to 200 for 
% autorotation and airspeed of 92km/h then gives a power-off descent 
% rate of 366 to 427m/min. 
% 
% Fatigue AGARD
% http://www.dtic.mil/dtic/tr/fulltext/u2/737398.pdf
% http://ftp.rta.nato.int/public//PubFullText/AGARD/AG/AGARD-AG-292///AGARD-AG-292.pdf
% http://www.wattflyer.com/forums/showthread.php?t=52194&page=22 B-47 thread
% 
% http://www.youtube.com/watch?v=dlYM7Y8Y8UU
% http://www.youtube.com/watch?v=mdP1jCWF3Eo
% http://www.rotorspot.nl/n-16.php
%
% Kellet wind tunnel testz
% http://crgis.ndc.nasa.gov/historic/643_Test_85_-_YG-1A_Autogiro
% http://rotored.arc.nasa.gov/
% https://www.youtube.com/watch?feature=player_embedded&v=KwwENgQ9reY
% 
% tradewind
% http://www.fspilotshop.com/virtavia-r3y-tradewind-for-fsx-p-4489.html?sid=e13221e774490586766116bab8165f05&utm_campaign=july12&utm_medium=email&utm_source=newsletter
% Hofstra
% http://www.hofstra.edu/community/culctr/autgir/autgir_photos.html
% Arliss Riggs
% http://www7a.biglobe.ne.jp/~GYROS/011022.html
% 
% Cierva Patent Collection 
% http://www.rotaryforum.com/forum/showthread.php?t=35602
%
% http://www.flightglobal.com/pdfarchive/view/1953/1953%20-%200094.html
% 
% One off Stinson
% http://www.flickr.com/photos/flyermedia/5141457126/
% http://www.flickr.com/photos/flyermedia/5140853521/
% http://www.flickr.com/photos/flyermedia/5140853099/in/photostream/
% R-4 fuselage
% http://www.flickr.com/photos/tom-margie/2130506565/sizes/z/in/set-72157603530016343/
% Wildcat
% http://www.flickr.com/photos/tom-margie/2133737753/sizes/z/in/set-72157603530016343/
% 
% Werkstattpraxis fuer den Bau von Gleit- und Segelflugzeugen
% http://www.lu5dje.com.ar/manual9.pdf
% http://www.holzleicht-flugzeugbau.de/Heuser_Dateien/H-IV/Prj_H-IV.html
% http://www.rotaryforum.com/forum/archive/index.php/t-22278.html
% Houston Modeling Of Helicopters
% http://www.google.de/url?sa=t&rct=j&q=0%20houston%2Cs.%20s.%2C%20%E2%80%9Cmodeling%20and%20analysis%20of%20helicopter%20flight%20mechanics%20in%20autorotation%2C&source=web&cd=11&cad=rja&ved=0CCkQFjAAOAo&url=http%3A%2F%2Fwww.tech-domain.com%2Fforum.php%3Fmod%3Dattachment%26aid%3DMTQ5NDZ8NjEwN2FjMmN8MTM3MTg3MzgxNHwwfDMzNzE1&ei=mKftUevND43AswbQhoHoBw&usg=AFQjCNEwwll_2rGCwfXBE8PzZE08zvGTIw&bvm=bv.49478099,d.ZWU
% https://www.youtube.com/watch?feature=player_embedded&v=KwwENgQ9reY
% FICON
% http://www.youtube.com/watch?v=D3bECsGocqE   1
% http://www.youtube.com/watch?v=IvEks8gbmgI   2
% B-36
% http://www.youtube.com/watch?v=7i92hwEDW58
% B-47
% http://www.youtube.com/watch?v=6cIgTAtj4E4
% P-26
% http://www.youtube.com/watch?v=3zXkVQnVmuo
% Jap Wings
% http://www.youtube.com/watch?v=v1iVfrQUQQM
% http://www.britishpathe.com/video/jet-helicopter/query/helicopter
% http://www.lowtechmagazine.com/2012/09/jobs-of-the-future-cargo-cyclist.html
% http://www.notechmagazine.com/2013/07/from-europe-to-america-by-sail.html
% http://www.classicbikersclub.com/articles/2013-04/yesterdays-antique-motorcycles
% http://www.driver.de/fiat-124-spider-ferrari-fuer-normalverbraucher/id_63021786/index
% Quattar MTO flying 25knots wind
% http://www.youtube.com/watch?v=2tXqWb63UdQ
% 
% Ice Pilots
% http://www.youtube.com/watch?feature=player_embedded&v=A6MjqDb-FT8&t=58
% Bush Pilots
% http://www.youtube.com/watch?v=vL3b-qzjo2U
% http://aerosociety.com/Assets/Docs/Publications/The%20Journal%20of%20Aeronautical%20History/2012-07_Evolution_of_Rotorcraft-Gibbings.pdf
%
% Guided Samara:
% http://dspace.mit.edu/bitstream/handle/1721.1/42047/229893867.pdf?sequence=1

