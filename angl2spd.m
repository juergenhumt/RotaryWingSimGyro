function [vH, vV, vHmph, vVfpm]=angl2spd(gamC, mu, alfaS, thtaF, a1, nRot,rRot)
   oM= pi*nRot/30;
   vD= mu*oM*rRot/cos(deg2rad(alfaS + a1));
   
   vH=vD*cos(thtaF+a1);
   vHmph=vH*3.6/1.605;
   vV=vH*tan(deg2rad(gamC));
   vVfps=vV/0.3048;
   vVfpm=vVfps*60;
   
   fprintf('vH %8.3f  vV %8.3f  vH %8.1f [mph]   vV %8.1f [fpm]\n',vH,vV,vHmph,vVfpm);
   % https://www.youtube.com/watch?v=hUsOQyfnfDY
   
   