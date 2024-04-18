%                                                          1    2    3    4    5     6        7       8      9
function [fR, mR, oMDotMR, fZMR] = AcFrcMmtSummationOutF( dlc, dle, dla, dlp, dlt, phiFus, thtaFus, oMTrim, Qmr, vErth, omB, fZMR,  outP, outFile)
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
%
% This routine calculates resultant forces and moments of the whole model 
% for trim calculation and perturbation values for selected coordinates
% to calculate stability derivatives. Routine is not part of NACA tm 73254
%
% Input:
%   1) dlc: collective input
%   2) dle: longitudinal cyclic input
%   3) dla: lateral cyclic input
%   4) dlp: pedal (tail rotor collective/ rudder) input
%   5) dlt: throttle setting
%   6) phiFus : trim roll angle
%   7) thtaFus: trim pitch angle 
%   8) oMR  : rotor tip speed
%   9) Qmr  : main rotor torque
%      vErth: velocities in earth coordinates
%      omC  : angular velocities body coordinates
%      outP : 0 < outP < 1000 -> print results
%        > 1000 return perturbation value for selected
%        coordinate to calculate stability derivatives
%
% Output
% fR      : resultant forces
% mTR     : resultant moments
% oMdotMR : rotor speed variation
% fZMR    : z component of rotor force
%
global rhoAir fpA A1STrim B1STrim A1CL B1CL hMR lVF C5 C6 C7 thta1 epsGlo iRigMR...
       mHeliG mHeli Iyy sigma gamRot rRot aLift thtaClmb iRigHS lmd xcg xCG...
       dDrgN dDrg1 dDrg2 f2m bTl rRot logFile psiTrim lbForce gloVars...
       aNCff a1Cff b1Cff laCff cTCff rrpmCff rhoAir aDiskMR dlpArry dlpSlct


 pB = omB(1);  qB=omB(2); rB = omB(3);
 zB = 1.0e4;

 vBdy = RTrfFunc(phiFus, thtaFus, psiTrim, vErth(1), vErth(2), vErth(3));
 
 if ((outP > 1000) & (outP < 2000))
   jD=outP-1000;
   if (2==jD)
     omB(jD)=omB(jD)+epsGlo;  
   else
     vBdy(jD)=vBdy(jD)+epsGlo;
   end  
 end
 thtaN = C5*dlc;
 AB1CP = d2ABCP(dle,dla);

 % control input angles and rigging angles are added
 % to give the control plane (swashplate) angles
 A1S= AB1CP(1) + A1CL;  
 B1S= AB1CP(2) + B1CL;

 % if outP > 10
   A1STrim = A1S;  % AB1CP   - 2.8
   B1STrim = B1S;  % AB1CP   + 4.6
 % end
  fRes=[0 0 0 0 0];

%  if abs(imag(thtaN)) > 0
%     fprintf('thtaN %10.5e\',thtaN);    
%  end
% 
  
%  if abs(imag(lamNf)) > 0
%     fprintf('lamNf %10.5e\',lamNf);    
%  end
% 
%  vBdyFus=vBdy;
%  vBdyFus(3)= vBdyFus(3) - lamI*oMTrim*rRot; % *cos(atan(vBdy(3)/vBdy(1)));
%

  
  [lamNf, xDotC, omDotC, anglC] = MainRotorCntrlFunc(vBdy, omB, oMTrim, thtaN, A1S, B1S, fZMR);
  lamI=anglC(6);

  [fFus, mFus] = FuselageFrcMmtFunc( vBdy, omB, anglC);
  [fHS, mHS] = HorizStabFrcMmtFunc(dle,  vBdy, omB, lamI, oMTrim, anglC);
  [fWN, mWN] = WingFrcMmtFunc(vBdy, omB, oMTrim, anglC); 
  [fPR, mPR] = PropellerFrcMmtFunc( vBdy, dlt);

  [fMR, mMR, oMDotMR] = MainRotorFrcMmtFunc(xDotC,omDotC, Qmr, anglC);
  fZMR=fMR(3);
  
  dlpArry(dlpSlct)=dlp;
  [fFN, mFN] = FinFrcMmtFunc( vBdy, omB, dlpArry(1));
  [fTR, mTR] = TailRotorFrcMmtFunc(dlpArry(2), vBdy, omB);
 
  % transformation of gravity to body coordinates
  frcG=frcGFunc(phiFus, thtaFus);

  %     1      2     3     4     5     6
  fR = fMR + fPR + fFus + fTR + fHS + fFN + frcG + fWN;
  mR = mMR + mPR + mFus + mHS + mFN + mTR + mWN;

 if (900==outP)
   options = optimset('MaxFunEvals',5000,'TolFun',1.0e-24,'Display','off');
   gloVars.fMRGlo = fMR;   gloVars.fFusGlo = fFus;   gloVars.fTRGlo = fTR;   gloVars.fHSGlo = fHS;   
   gloVars.fFNGlo = fFN;   gloVars.frcGGlo = frcG;   gloVars.fPRGlo = fPR;
   gloVars.mMRGlo = mMR;   gloVars.mFusGlo = mFus;   gloVars.mHSGlo = mHS;   gloVars.mFNGlo = mFN;   
   gloVars.mTRGlo = mTR;   gloVars.mPRGlo  = mPR; gloVars.torque=Qmr;
   gloVars.QR     = oMDotMR; gloVars.thtaN=anglC(5); gloVars.mu=anglC(7);
   [xRes fVal] = fsolve(@oMitr,oMTrim, options, vBdy, omB, thtaN, A1S, B1S, -fZMR);
   gloVars.DoM_it= (xRes - oMTrim); gloVars.alfaNf= anglC(13); gloVars.a1=anglC(1); gloVars.cTs=anglC(14)/sigma;
 end

 if (outP > 1000)
   options = optimset('MaxFunEvals',5000,'TolFun',1.0e-24,'Display','off');   
   outS='xxxx';
   if 1==jD
     outS='u - derivatives';
     [xRes fVal] = fsolve(@oMitr,oMTrim, options, vBdy, omB, thtaN, A1S, B1S, -fZMR);
     gloVars.DoMDu_it= (xRes - oMTrim)/epsGlo;
   end
   if 2==jD
     outS='q - derivatives';
     [xRes fVal] = fsolve(@oMitr,oMTrim, options, vBdy, omB, thtaN, A1S, B1S, -fZMR);
     gloVars.DoMDq_it= (xRes - oMTrim)/epsGlo;    
   end
   if 3==jD
     outS='w - derivatives';
     [xRes fVal] = fsolve(@oMitr,oMTrim, options, vBdy, omB, thtaN, A1S, B1S, -fZMR);
     gloVars.DoMDw_it= (xRes - oMTrim)/epsGlo;    
   end
   
   fprintf(outFile,'\n ====== %s ===== \n',outS);

   
   
   fMR= (fMR -  gloVars.fMRGlo)/mHeli/epsGlo;   fFus= (fFus -  gloVars.fFusGlo)/mHeli/epsGlo;   fTR= (fTR -  gloVars.fTRGlo)/mHeli/epsGlo;   fHS= (fHS -  gloVars.fHSGlo)/mHeli/epsGlo;   
   fFN= (fFN -  gloVars.fFNGlo)/mHeli/epsGlo;   frcG= (frcG -  gloVars.frcGGlo)/mHeli/epsGlo; fPR= (fPR -  gloVars.fPRGlo)/mHeli/epsGlo;
   mMR= (mMR-  gloVars.mMRGlo)/Iyy/epsGlo; mFus= (mFus-  gloVars.mFusGlo)/Iyy/epsGlo;   mTR= (mTR-  gloVars.mTRGlo)/Iyy/epsGlo;   mHS= (mHS-  gloVars.mHSGlo)/Iyy/epsGlo;   
   mFN= (mFN-  gloVars.mFNGlo)/Iyy/epsGlo; mPR=(mPR-gloVars.mPRGlo)/Iyy/epsGlo;
   oMDotMR = (oMDotMR - gloVars.QR);

   fprintf(outFile,'fMR(1)  %15.6f   fMR(2)  %15.6f   fMR(3)  %15.6f\n',fMR(1),fMR(2),fMR(3));
   fprintf(outFile,'fTR(1)  %15.6f   fTR(2)  %15.6f   fTR(3)  %15.6f\n',fTR(1),fTR(2),fTR(3)); 
   fprintf(outFile,'fFus(1) %15.6f   fFus(2) %15.6f   fFus(3) %15.6f\n',fFus(1),fFus(2),fFus(3));
   fprintf(outFile,'fFN(1)  %15.6f   fFN(2)  %15.6f   fFN(3)  %15.6f\n',fFN(1),fFN(2),fFN(3));
   fprintf(outFile,'fHS(1)  %15.6f   fHS(2)  %15.6f   fHS(3)  %15.6f\n',fHS(1),fHS(2),fHS(3));
   fprintf(outFile,'fPR(1)  %15.6f   fPR(2)  %15.6f   fPR(3)  %15.6f\n',fPR(1),fPR(2),fPR(3));
   fprintf(outFile,'fCG(1)  %15.6f   fCG(2)  %15.6f   fCG(3)  %15.6f\n',frcG(1),frcG(2),frcG(3));
   
   fprintf(outFile,'fR(1)  %15.7e   fR(2)  %15.7e   fR(3)  %15.7e\n',fR(1),fR(2),fR(3));
  
   fprintf(outFile,'\n Moment Equilibrium %7.0f \n',outP);
   fprintf(outFile,'mMR(1)  %15.6f   mMR(2)  %15.6f   mMR(3)  %15.6f\n',mMR(1),mMR(2),mMR(3));
   fprintf(outFile,'mTR(1)  %15.6f   mTR(2)  %15.6f   mTR(3)  %15.6f\n',mTR(1),mTR(2),mTR(3));
   fprintf(outFile,'mFus(1) %15.6f   mFus(2) %15.6f   mFus(3) %15.6f\n',mFus(1),mFus(2),mFus(3));
   fprintf(outFile,'mFN(1)  %15.6f   mFN(2)  %15.6f   mFN(3)  %15.6f\n',mFN(1),mFN(2),mFN(3));
   fprintf(outFile,'mHS(1)  %15.6f   mHS(2)  %15.6f   mHS(3)  %15.6f\n',mHS(1),mHS(2),mHS(3));
   fprintf(outFile,'mPR(1)  %15.6f   mPR(2)  %15.6f   mPR(3)  %15.6f\n',mPR(1),mPR(2),mPR(3));

   fprintf(outFile,'mR(1)  %15.7e   mR(2)  %15.7e   mR(3)  %15.7e\n',mR(1),mR(2),mR(3));

   fprintf(outFile,'QR %15.7e\n',oMDotMR);

   xh= mMR(2) + mTR(2) + mFus(2)+ mFN(2) + mHS(2) + mPR(2);
   if 1==jD
     fprintf(outFile,'mU %15.7e\n',xh);
   end
   
   if 3==jD
     fprintf(outFile,'mW %15.7e\n',xh);
   end

   
 end
 
 if ((outP > 0) & (outP < 999))
%   [cTs, cT, lamNf, lamI, alfaD, alfaDeg,alfaNf]  = clcRotState2(thtaN, norm(vBdy)/(oMTrim*rRot));

   
   fprintf('\n==a=a=== frc func ===a=a==\n');
   a1  = anglC(1);
   b1  = anglC(3);

   lamNf = anglC(12);
   lamI  = anglC(6);
   
   thtaN = anglC(5);
   
   cT = anglC(15);
   cTs = cT/sigma;
   
   mu=anglC(7);    mu2=mu*mu;  thta75=anglC(5) + 0.75*thta1;
   
   oM  = anglC(8);
   oMR = oM*rRot;
   
   aNBr = 0.125*gamRot*( thta75*(1+mu2) + 4/3*lamNf);
   a1Br = 2*mu*(4/3*thta75+lamNf)/(1-0.5*mu2)*(1+1.5*mu);
   
   % anglC(6) is induced inflow
   b1Br= 4/3*(mu*aNBr + 1.1*anglC(6))/(1+0.5*mu2);

   lamD = lamNf + mu*a1; 
   alfaNfAux = rad2deg(anglC(13));
   
   cTHov=mHeliG/(rhoAir*pi*rRot^2*(oMTrim*rRot)^2);
   laHov=sqrt(0.5*cTHov);
   
   alfD = rad2deg(atan((lamD + anglC(6))/mu));
   alfaNf =alfD- rad2deg(a1);
   cL =2*cT*cos(deg2rad(alfaNf))^3/mu^2;
   cLs = cL/sigma;
   xh = 2*cT/sigma/aLift;
   cT2sa = eq6(mu,lamNf,thtaN);
   xh2 = cT2sa/mu^2*aLift;
   cTs2 = 0.5*aLift*cT2sa;
   
   aN475 = crvVal(mu, aNCff);
   a1475 = crvVal(mu, a1Cff);
   b1475 = crvVal(mu, b1Cff);
   la600 = crvVal(mu, laCff);
   rrpm600 = crvVal(mu,rrpmCff);
   xh = mu*pi*rrpm600/30*rRot;
   cT600 = crvVal(xh,cTCff);
   
   
   alfaNfaux= la600/mu + cT600/(2*mu*sqrt(mu^2 + la600^2));
   
   fprintf('xMR %8.3f    xcg %8.3f    dxCG %8.3f    iRigMR %6.1f\n',lmd/f2m,xcg/f2m,xCG/f2m,rad2deg(iRigMR));
   fprintf('cTs %7.5f   cT %7.5f   cLs %7.5f   cL %7.5f\n',cTs,cT,cL/sigma,cL);
   cTs_cr=0.29*thtaN; 
   fprintf('mu %6.3f  mu/laHov %7.5f  thtaN/cTs %7.5f  cTs_cr %7.5f  cTs_cr/cTs %7.5f\n',mu,mu/laHov,thtaN/cTs,cTs_cr,cTs_cr/cTs);
   
   qRot= gloVars.torque;
   P_kW=qRot*oMTrim/1000;     P_hp=P_kW/0.743;
   
   fprintf('Qrot[Nm]  %9.1f   P[kW] %7.2f   P[hp] %7.2f          dQ %9.1e\n',qRot,P_kW,P_hp,anglC(14));
   fprintf('aN  %10.5f°    a1 %10.5f°   b1 %10.5f°\n',aNBr,a1Br,b1Br);
   fprintf('Br   aN %8.2f°       a1 %7.2f°      b1 %7.2f°\n',rad2deg(aNBr),rad2deg(a1Br),rad2deg(b1Br));
   fprintf('716  aN %8.2f°       a1 %7.2f°      b1 %7.2f°\n',rad2deg(anglC(9)),rad2deg(anglC(1)),rad2deg(anglC(3)));
   a1475 = rad2deg(eq2(mu,-0.056,thtaN*9.4/7));
   fprintf('mmt  aN %8.2f°       a1 %7.2f°      b1 %7.2f°\n',aN475,a1475,b1475);
   fprintf('A1S %10.5f°    B1S %10.5f°\n',anglC(2),anglC(4));
   fprintf('A1S %7.2f°       B1S %7.2f°\n',rad2deg(anglC(2)),rad2deg(anglC(4)));
   
   fprintf('oM %6.3f[1/s]   oMR %6.1f[m/s]   oMR %6.1f[ft/s]    n %5.1f   nmmt %5.1f [1/min]\n',oMTrim,oMTrim*rRot,oMTrim*rRot/f2m,oMTrim/(2*pi)*60,crvVal(norm(vErth)*3.6/1.605,rrpmCff))
   fprintf('lamD %9.5f     lamI %9.5f     lamNf %10.5f     lamNf/laHov  %10.5f     lammt %10.5f\n',lamD,anglC(6),lamNf,lamNf/laHov,la600);
   xh = iRigMR - B1S + a1 + thtaFus;

   a1s = - B1S + a1;
   %       = alfaNf    + B1
   alfaS_1 = (anglC(13) + B1S);
   alfaS_2 = (-thtaClmb + thtaFus + iRigMR);
   
   xh1 = rad2deg((mu*(8/3*thtaN + 2*(thta1 + mu*alfaS_1 - anglC(6)*oMTrim*rRot)))/(1-0.5*mu*mu));
   xh2 = rad2deg((1+1.5*mu*mu)*B1S/(1-0.5*mu*mu));
   %                                               alfaI
   a1s_th = (mu*(8/3*thtaN + 2*(thta1 + mu*alfaS_1 - anglC(4)*oMTrim*rRot)) - (1+1.5*mu*mu)*B1S)/(1-0.5*mu*mu);
   
   alfaS_1 = rad2Deg(alfaS_1);
   alfaS_2 = rad2Deg(alfaS_2);
   
   alfaD=anglC(13)+anglC(1);
   
   
   fprintf('alfaNf %5.2f°  alfaD   %5.2f°   alfaS %5.2f°   als %8.2f°   als_th %8.2f°   alNfmmt %8.2f° \n',rad2deg(anglC(13)),rad2deg(alfaD),...
       rad2deg(anglC(13)+anglC(1))-rad2deg(anglC(1)-anglC(4)),rad2deg(a1s),rad2deg(a1s_th),rad2deg(atan(alfaNfaux)));

   muZ = -xDotC(3)/(oMTrim*rRot);
   xCTNN=rhoAir*aDiskMR*sigma*(oMTrim*rRot)^2;
   cTs= -fZMR/xCTNN;
   
   thtaTR=C6*dlp+C7;
   fprintf('thtaN %6.2f°  thta1 %6.2f°  thtaClmb %6.2f°     thtaTR %8.2f°\n',rad2deg(thtaN),rad2deg(thta1),rad2deg(thtaClmb),rad2deg(thtaTR));
   fprintf('Trim   phi %8.2f°    thta %8.2f°    psi %8.2f°\n',rad2deg(phiFus),rad2deg(thtaFus),rad2deg(psiTrim));
   fprintf('dDrgN  %7.4f  dDrg1  %7.4f  dDrg2  %7.4f\n',dDrgN,dDrg1,dDrg2);
   xh=norm(vBdy)/f2m;
   fprintf('vBx %8.2f    vBy %8.2f    vBz %8.2f [m/s]   V %8.2f [fps]\n',vBdy(1),vBdy(2),vBdy(3),xh);
   fprintf('vEx %8.2f    vEy %8.2f    vEz %8.2f [fpm]   V %8.2f [mph]\n',vErth(1)*3.6/1.605,vErth(2),-60*vErth(3)/f2m, norm(vErth)*3.6/1.605);

   fhiMR=atan(fMR(1)/fMR(3));
   xh = rad2deg(fhiMR);
   dxMR= xCG - hMR*tan(fhiMR);
   xh2 = rad2deg(atan(gloVars.Fx/gloVars.Fz));
   % dxMR < 0 --> CoG behind      rotor force
   % dxMR > 0 --> CoG in front of rotor force
%    PvL=0.27;
%    [alfaPr, cTs, DvLN] = PvLMuIt(mu, thtaN, vErth, PvL,oMTrim);
%    cQ  = eq9a(mu,lamNf,thtaN,oM) - eq11d(mu,lamNf,thtaN,oMTrim);
%    PvL = cQ/mu/(cTs2*sigma);
%    fprintf('cT %7.4f   cQ %8.5f PvL %8.5f \n',cTs2*sigma,cQ,PvL);
%    return 
%    
%    
heliDeriv3104(alfaD, thtaN, mu, muZ, oMTrim, cTs, lamD, lamI, a1, b1, a1s, vBdy, anglC, outFile);   
     
cT2sa = eq6(mu, lamNf, thtaN);
cQ2s = eq11d(mu, lamNf, thtaN) - eq9a(mu, lamNf, thtaN);
DvLp = fpA*rhoAir/2*norm(vErth)^2/mHeliG;
% 0.275
PvL  = cQ2s/(mu*cT2sa*aLift);
DvLN = abs(eq13(mu,lamNf,thtaN)/(mu*cT2sa));
fprintf('PvL %4.2f   DvLN %9.5f    lamNf %9.5f\n',PvL,DvLN,lamNf);
     
     
     
%   end 
   fprintf('MR     Fx %8.2f   Fz %8.2f   Fx/Fz %6.1f°    Fx/Fz2 %6.1f°    dxMR  %8.4f\n\n',gloVars.Fx,gloVars.Fz,xh,xh2,dxMR);
   % fMR2 = RTrfY(-iRigMR,gloVars.Fx,0,gloVars.Fz);
   fprintf('MR dsk Fx %8.2f    Fz %8.2f\n',gloVars.FxDsk,gloVars.FzDsk);
   
   fprintf('fMR(1)  %15.2f   fMR(2)  %15.2f   fMR(3)  %15.2f\n',fMR(1),fMR(2),fMR(3));
   fprintf('fTR(1)  %15.2f   fTR(2)  %15.2f   fTR(3)  %15.2f\n',fTR(1),fTR(2),fTR(3)); 
   fprintf('fPR(1)  %15.2f   fPR(2)  %15.2f   fPR(3)  %15.2f\n',fPR(1),fPR(2),fPR(3)); 
   fprintf('fFus(1) %15.2f   fFus(2) %15.2f   fFus(3) %15.2f\n',fFus(1),fFus(2),fFus(3));
   fprintf('fFN(1)  %15.2f   fFN(2)  %15.2f   fFN(3)  %15.2f\n',fFN(1),fFN(2),fFN(3));
   fprintf('fHS(1)  %15.2f   fHS(2)  %15.2f   fHS(3)  %15.2f\n',fHS(1),fHS(2),fHS(3));
   fprintf('fCG(1)  %15.2f   fCG(2)  %15.2f   fCG(3)  %15.2f\n',frcG(1),frcG(2),frcG(3));
   fprintf('fWN(1)  %15.2f   fWN(2)  %15.2f   fWN(3)  %15.2f\n',fWN(1),fWN(2),fWN(3));
   
   fprintf('fR(1)  %15.7e   fR(2)  %15.7e   fR(3)  %15.7e\n',fR(1),fR(2),fR(3));
   fprintf('fR(3)  %15.1f   fW(3)  %15.1f   [lb]',fMR(3)/lbForce,fWN(3)/lbForce);
     
   fprintf('\nMoment Equilibrium\n')
   fprintf('mMR(1)  %15.2f   mMR(2)  %15.2f   mMR(3)  %15.2f\n',mMR(1),mMR(2),mMR(3));
   fprintf('mTR(1)  %15.2f   mTR(2)  %15.2f   mTR(3)  %15.2f\n',mTR(1),mTR(2),mTR(3));
   fprintf('mPR(1)  %15.2f   mPR(2)  %15.2f   mPR(3)  %15.2f\n',mPR(1),mPR(2),mPR(3));
   fprintf('mFus(1) %15.2f   mFus(2) %15.2f   mFus(3) %15.2f\n',mFus(1),mFus(2),mFus(3));
   fprintf('mFN(1)  %15.2f   mFN(2)  %15.2f   mFN(3)  %15.2f\n',mFN(1),mFN(2),mFN(3));
   fprintf('mHS(1)  %15.2f   mHS(2)  %15.2f   mHS(3)  %15.2f\n',mHS(1),mHS(2),mHS(3));
   fprintf('mWN(1)  %15.2f   mWN(2)  %15.2f   mWN(3)  %15.2f\n',mWN(1),mWN(2),mWN(3));

   fprintf('mR(1)  %15.7e   mR(2)  %15.7e   mR(3)  %15.7e\n',mR(1),mR(2),mR(3));
   fprintf('==x=x=== end frc func ===x=x==\n');

end



if (outP == 999)
    
    
   [cTs, cT, lamNf, lamI, alfaD, alfaDeg,alfaNf]  = clcRotState2(thtaN, norm(vBdy)/(oMTrim*rRot));
    
   fprintf(outFile,'\n\n=I=I=I=I=I= Iteration Results =I=I=I=I=I=\n');
   mu=anglC(7);    mu2=mu*mu;  thta75=anglC(5) + 0.75*thta1;
   aNBr=0.125*gamRot*( thta75*(1+mu2) + 4/3*lamNf);
   a1Br= 2*mu*(4/3*thta75+lamNf)/(1-0.5*mu2)*(1+1.5*mu);
   a1  = anglC(1);
   
   % anglC(6) is induced inflow
   b1Br= 4/3*(mu*aNBr + 1.1*anglC(6))/(1+0.5*mu2);
   b1  = anglC(3);
   
   % cTs = anglC(9);
   lamD=lamNf+mu*a1; 
   % xh=sigma*anglC(9)/(2*bTl^2*sqrt(mu2+lamD*lamD));
   % alfD= rad2deg(atan((lamD + xh)/mu));
   % cT = cTs*sigma;
   
   alfD = rad2deg(atan((lamD + anglC(6))/mu));
   
   alfaNf=alfD- rad2deg(a1);
   cL=2*cT*cos(deg2rad(alfaNf))^3/mu^2;
   fprintf(outFile,'xMR %8.3f    xcg %8.3f    dxCG %8.3f    iRigMR %6.1f\n',lmd/f2m,xcg/f2m,xCG/f2m,rad2deg(iRigMR));
   fprintf(outFile,'cTs %7.5f   cT %7.5f   cLs %7.5f   cL %7.5f\n',cTs,cT,cL/sigma,cL);
   fprintf(outFile,'aN  %10.5f°    a1 %10.5f°   b1 %10.5f°\n',aNBr,a1Br,b1Br);
   fprintf(outFile,'Br   aN %8.2f°       a1 %7.2f°      b1 %7.2f°\n',rad2deg(aNBr),rad2deg(a1Br),rad2deg(b1Br));
   fprintf(outFile,'716  aN %8.2f°       a1 %7.2f°      b1 %7.2f°\n',rad2deg(aNBr),rad2deg(anglC(1)),rad2deg(anglC(3)));

   fprintf(outFile,'mu %9.6f lamD %10.5f    lamI %10.5f     mu*a1 %10.5f     lamNf %10.5f\n',mu,lamD,anglC(6),mu*a1,lamNf);
   xh = iRigMR - B1S + a1 + thtaFus;

   xh = - B1S + a1;
   
   fprintf(outFile,'A1S %10.5f°    B1S %10.5f°\n',anglC(2),anglC(4));
   fprintf(outFile,'A1S %8.2f°    B1S %8.2f°    alfD %8.2f°   alfaNf %8.2f°   alSm %8.2f°\n',rad2deg(anglC(2)),rad2deg(anglC(4)),alfD,alfaNf,rad2deg(xh));
   thtaTR=C6*dlp+C7;
   fprintf(outFile,'thtaN  %8.2f°    thtaTR %8.2f°\n',rad2deg(thtaN),rad2deg(thtaTR));
   fprintf(outFile,'Trim   phi %8.2f°    thta %8.2f°\n',rad2deg(phiFus),rad2deg(thtaFus));
   fprintf(outFile,'vBx %8.2f    vBy %8.2f    vBz %8.2f\n',vBdy(1),vBdy(2),vBdy(3));
  
   fprintf(outFile,'fMR(1)  %15.2f   fMR(2)  %15.2f   fMR(3)  %15.2f\n',fMR(1),fMR(2),fMR(3));
   fprintf(outFile,'fTR(1)  %15.2f   fTR(2)  %15.2f   fTR(3)  %15.2f\n',fTR(1),fTR(2),fTR(3)); 
   fprintf(outFile,'fPR(1)  %15.2f   fPR(2)  %15.2f   fPR(3)  %15.2f\n',fPR(1),fPR(2),fPR(3)); 
   fprintf(outFile,'fFus(1) %15.2f   fFus(2) %15.2f   fFus(3) %15.2f\n',fFus(1),fFus(2),fFus(3));
   fprintf(outFile,'fFN(1)  %15.2f   fFN(2)  %15.2f   fFN(3)  %15.2f\n',fFN(1),fFN(2),fFN(3));
   fprintf(outFile,'fHS(1)  %15.2f   fHS(2)  %15.2f   fHS(3)  %15.2f\n',fHS(1),fHS(2),fHS(3));
   fprintf(outFile,'fCG(1)  %15.2f   fCG(2)  %15.2f   fCG(3)  %15.2f\n',frcG(1),frcG(2),frcG(3));
   fprintf(outFile,'fWN(1)  %15.2f   fWN(2)  %15.2f   fWN(3)  %15.2f\n',fWN(1),fWN(2),fWN(3));
   
   fprintf(outFile,'fR(1)  %15.7e   fR(2)  %15.7e   fR(3)  %15.7e\n',fR(1),fR(2),fR(3));
  
   fprintf(outFile,'\nMoment Equilibrium\n');
   fprintf(outFile,'mMR(1)  %15.2f   mMR(2)  %15.2f   mMR(3)  %15.2f\n',mMR(1),mMR(2),mMR(3));
   fprintf(outFile,'mTR(1)  %15.2f   mTR(2)  %15.2f   mTR(3)  %15.2f\n',mTR(1),mTR(2),mTR(3));
   fprintf(outFile,'mPR(1)  %15.2f   mPR(2)  %15.2f   mPR(3)  %15.2f\n',mPR(1),mPR(2),mPR(3));
   fprintf(outFile,'mFus(1) %15.2f   mFus(2) %15.2f   mFus(3) %15.2f\n',mFus(1),mFus(2),mFus(3));
   fprintf(outFile,'mFN(1)  %15.2f   mFN(2)  %15.2f   mFN(3)  %15.2f\n',mFN(1),mFN(2),mFN(3));
   fprintf(outFile,'mHS(1)  %15.2f   mHS(2)  %15.2f   mHS(3)  %15.2f\n',mHS(1),mHS(2),mHS(3));
   fprintf(outFile,'mWN(1)  %15.2f   mWN(2)  %15.2f   mWN(3)  %15.2f\n',mWN(1),mWN(2),mWN(3));

   fprintf(outFile,'mR(1)  %15.7e   mR(2)  %15.7e   mR(3)  %15.7e\n',mR(1),mR(2),mR(3));
   fprintf(outFile,'==x=x===x===x=x==\n');
   fhiMR=atan(fMR(1)/fMR(3));
   xh = rad2deg(fhiMR);
   dxMR= xCG - hMR*tan(fhiMR);

   % dxMR < 0 --> CoG  behind     rotor force
   % dxMR > 0 --> CoG in front of rotor force
   fprintf(outFile,'MR Fx %8.2f    Fz %8.2f   thtaClmb %8.2f   torque %8.2f\n',gloVars.Fx,gloVars.Fz,rad2deg(thtaClmb),gloVars.torque);
   fprintf(outFile,'dxMR  %8.4f  Fx/Fz %6.1f°        mu %9.5f   oM %5.1f  n %6.1f\n',dxMR,xh,mu,oMTrim,oMTrim/(2*pi)*60);   
   fprintf(outFile,'\n=I=I=I=I=I= End Of Iteration Results =I=I=I=I=I=\n\n');

end


