function trimOut = acTrim3DFM(trmVr, vErth, omBd, jOut, cDmp, outFile);
%
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
% This routine trims the aircraft to zero translational and rotational accelerations 
% and zero rotational velocities. To this end small perturbations are imposed on the 
% trim variables and the jacobian that results is inverted to calculate an improved 
% value. It is basicaly a Newton search algorithm. The routine is based on the one 
% proposed by Mark E. Dreier in "Introduction to Helicopter and Tiltrotor Flight
% Dynamics"
%
% 
% Input:
% trimVr: the eight trim variables. These are
% 1)  dlc: main rotor collective 
% 2)  dle: longitudinal stick input      
% 3)  dla: lateral stick input
% 4)  dlp: tail rotor collective via pedals
% 5)  dlThrtl : throttle setting
% 6)  thtaTrim: pitch angle
% 7)  phiTrim : roll angle
% 8)  oM      : rotor speed
% vErth : translational velocities in earth coordinates
% omBd: angular velocities in body coordinates
% oM  : rotationl speed of rotor
%
% Output
% trmVr: the varible values for trimmed flight 
% 
global mHeli mHeliG rRot Ixx Iyy Izz epsGlo thtaClmb...
  firstOutC firstOutB firstOutMR  firstOutW itrType rfrncVal gloVars

outP=0;

% trmVec=[1 2 3 4 5 6 7 8 9];
if (0==itrType)
% using trim variables 1 through 7 means that
% collective varies and rotor speed is fixed 
  fprintf('autogyro mode: iterating thtaN\n')
  trmVec=[1 2 3 4 5 6 7];
elseif (1==itrType)
% using trim variables 2 through 8 means that
% collective is fixed and rotor speed varies
  fprintf('autogyro mode: iterating oM\n')
  trmVec=[2 3 4 5 6 7 8];
  % trmVec=[2 3 5 6 7 8];
elseif (2==itrType)
% using trim variables 2-4/6-9 means that
% collective and rotor torque vary while 
% rotor speed is fixed, engine power is zero
  fprintf('helicopter mode: iterating thtaN\n')
  trmVec=[1 2 3 4 6 7 9];  
else
  fprintf('ERROR: Invalid itrType\nexiting from trim\n');
  return;
end 

 jTDim=length(trmVec);
 
 vN=zeros([1 3]);
 vNj=zeros([1 jTDim]);

 fRes=vN;
 mRes=vN;


 cont    = vNj;
 deltt6  = vNj;

 acc0Trim = vNj;
 accpTrim = vNj;
 accmTrim = vNj;
 deltt    = vNj;
 jcb3D    = zeros(jTDim,jTDim);
 jcb3Dinv = zeros(jTDim,jTDim);

%                        AcFrcMmtSummationOutF(dlc,      dle,      dla,      dlp,      dlt,     phi      thtaFus     oM       
 [fext mext oMDot fZMR]= AcFrcMmtSummationOutF(trmVr(1), trmVr(2), trmVr(3), trmVr(4), trmVr(5),trmVr(6), trmVr(7), trmVr(8), trmVr(9), vErth, omBd,  mHeliG, outP, outFile);
%      phiFus,    thtaFus, psiFus    oMTrim,  vErth, omB,   fZMR,   outP, outFile)

 fZMR=-fZMR;
 
 acc0Trim = [fext mext oMDot];
 anorm=norm(fext) + norm(mext);
	
 eps = 1.0e-6;   itEnd=  150;
	% #### Trim iter start ####
 iter= 1; epsEnd=1.0e-9;
 while ((anorm > epsEnd) & (iter < itEnd))
	% Now build the Jacobian. Perturb the elements of the trim vector trmVr.
	%
		for jac=1:jTDim 
		%
		% Fwd perturbation
		 trmVr(trmVec(jac)) = trmVr(trmVec(jac)) + eps;
		
		 [fext mext oMDot xh]= AcFrcMmtSummationOutF(trmVr(1), trmVr(2), trmVr(3),trmVr(4), trmVr(5), trmVr(6),  trmVr(7), trmVr(8), trmVr(9),vErth, omBd, fZMR, outP, outFile);
		  accpTrim = [fext mext oMDot];
		%
		% Bkwd perturbation
		 trmVr(trmVec(jac)) = trmVr(trmVec(jac)) - 2 * eps;
		 
		 [fext mext oMDot xh]= AcFrcMmtSummationOutF(trmVr(1), trmVr(2), trmVr(3),trmVr(4), trmVr(5), trmVr(6),  trmVr(7), trmVr(8), trmVr(9), vErth, omBd,  fZMR, outP, outFile);
		  accmTrim = [fext mext oMDot];
		
		%Return trim column to starting value in this iteration
		  trmVr(trmVec(jac)) = trmVr(trmVec(jac)) + eps;
		
        % Build the derivative column.
          for iac = 1:jTDim
             jcb3D(iac, jac) = 0.5*(accpTrim(iac) - accmTrim(iac))/eps;
          end   % iSlct(iac)
          
      end   % jac
	
      printJcb(jcb3D,jTDim,firstOutB);
	  jcb3Dinv = inv(jcb3D);
	  deltt= acc0Trim*jcb3Dinv';%'
 	  printJcb(jcb3Dinv,jTDim,firstOutB);
	  firstOutB=0;
	  
	  for i = 1:jTDim
        dlt=real(deltt(i));
        % dlt=deltt(i);
	    trmVr(trmVec(i)) = trmVr(trmVec(i)) - cDmp*dlt;
	  end   % i
	 
	%
	% Evaluate error 
	  [fext mext oMDot fZMR]= AcFrcMmtSummationOutF(trmVr(1), trmVr(2), trmVr(3),trmVr(4), trmVr(5), trmVr(6),  trmVr(7),  trmVr(8), trmVr(9), vErth, omBd,  fZMR, outP, outFile);
      acc0Trim = [fext mext oMDot];
      anorm=norm(fext) + norm(mext) + abs(oMDot);
      fZMR=-fZMR;

	  iter = iter + 1;  
  end % while

  if iter >= itEnd
     fprintf('+++ Error in actrim3D! +++\ntrim iter %d =itEnd    anorm %10.5e\n',iter,anorm);
     trimOut=zeros([1 8]);
     % return    
  else
     fprintf(outFile,'trim iter %d     anorm %10.5e\n',iter,anorm);
  end
  
  trimOut   = trmVr;
  trimOut(9)= fZMR;
  
  if jOut < 2
  outP= 900; firstOutC=1.0;  firstOutMR=1.0; firstOutW=1.0;
  [fextN mextN QRN]= AcFrcMmtSummationOutF(trmVr(1), trmVr(2), trmVr(3),trmVr(4), trmVr(5), trmVr(6),  trmVr(7),  trmVr(8), trmVr(9), vErth, omBd,  fZMR, outP, outFile);
  [cTs, cT, lamNf, vi, alfaD, alfaDeg,alfaNf]  = clcRotState2(gloVars.thtaN, gloVars.mu);
  firstOutC=0.0;
  
  outP=1001; % u + eps
  [fext1 mext1 QR1 TMR]= AcFrcMmtSummationOutF(trmVr(1), trmVr(2), trmVr(3),trmVr(4), trmVr(5), trmVr(6),  trmVr(7),  trmVr(8), trmVr(9), vErth, omBd,  fZMR, outP, outFile);
  
  outP=1002;  % q + eps
  [fext2 mext2  QR2]= AcFrcMmtSummationOutF(trmVr(1), trmVr(2), trmVr(3),trmVr(4), trmVr(5), trmVr(6),  trmVr(7),  trmVr(8), trmVr(9), vErth, omBd,  fZMR, outP, outFile);
  
  outP=1003;  % w + eps
  [fext3 mext3 QR3]= AcFrcMmtSummationOutF(trmVr(1), trmVr(2), trmVr(3),trmVr(4), trmVr(5), trmVr(6),  trmVr(7),  trmVr(8), trmVr(9), vErth, omBd,  fZMR, outP, outFile);

  outP= 999;  % oM + eps
  [fext8 mext8 QR8]= AcFrcMmtSummationOutF(trmVr(1), trmVr(2), trmVr(3),trmVr(4), trmVr(5), trmVr(6),  trmVr(7),  trmVr(8)+epsGlo, trmVr(9), vErth, omBd,  fZMR, outP, outFile);
  

  fprintf(outFile,'\n         ### stability derivatives ###\n');
  fprintf(outFile,' R ._. R ._. << RotaryWingSim Data >> ._. R ._. R\n\n');

  fprintf(outFile,'       u              w            q         oM\n');
  xU= (fext1(1)-fextN(1))/mHeli/epsGlo;    xW= (fext3(1)-fextN(1))/mHeli/epsGlo;  xQ= (fext2(1)-fextN(1))/mHeli/epsGlo;  xoM=(fext8(1)-fextN(1))/mHeli/epsGlo;
  zU= (fext1(3)-fextN(3))/mHeli/epsGlo;    zW= (fext3(3)-fextN(3))/mHeli/epsGlo;  zQ=(fext2(3)-fextN(3))/mHeli/epsGlo;  zoM=(fext8(3)-fextN(3))/mHeli/epsGlo;
  % mU=(mext1(2)-mextN(2))/Iyy/rRot/epsGlo; mW=(mext3(2)-mextN(2))/Iyy/rRot/epsGlo;    mQ=(mext2(2)-mextN(2))/Iyy/rRot/epsGlo;    moM=(mext8(2)-mextN(2))/Iyy/rRot/epsGlo;  
  mU= (mext1(2)-mextN(2))/Iyy/epsGlo; mW=(mext3(2)-mextN(2))/Iyy/epsGlo;    mQ=(mext2(2)-mextN(2))/Iyy/epsGlo;    moM=(mext8(2)-mextN(2))/Iyy/epsGlo;  
  
  fprintf(outFile,'x  %10.6f    %10.6f   %10.6f   %10.6f\n',xU,xW,xQ,xoM);
  fprintf(outFile,'z  %10.6f    %10.6f   %10.6f   %10.6f\n',zU,zW,zQ,zoM);
  fprintf(outFile,'m  %10.6f    %10.6f   %10.6f   %10.6f\n',mU,mW,mQ,moM);
  
    xh=(fext8(1)-fextN(1))/mHeli/epsGlo;
    gloVars.DxDoM= xh;
    
    zh=(fext8(3)-fextN(3))/mHeli/epsGlo;
    gloVars.DzDoM= zh;

    mh=(mext8(2)-mextN(2))/mHeli/epsGlo;
    gloVars.DmDoM=mh;

  fprintf(outFile,'\n      measured  data\n');
  fprintf(outFile,'       u                w              q\n');
  fprintf(outFile,'x  %10.6f     %10.6f      %10.6f\n',rfrncVal(1),rfrncVal(2),rfrncVal(3));
  fprintf(outFile,'z  %10.6f     %10.6f      %10.6f\n',rfrncVal(4),rfrncVal(5),rfrncVal(6));
  fprintf(outFile,'m  %10.6f     %10.6f      %10.6f\n\n',rfrncVal(7),rfrncVal(8),rfrncVal(9));


  QRU=QR1-QRN; QRW=QR2-QRN; QRQ=QR3-QRN; QRoM=QR8-QRN;
 
  fprintf(outFile,'Q  %10.6f     %10.6f      %10.6f      %10.6f\n\n',QRU,QRW,QRQ,QRoM);
  

  
if (jOut > 0.5) & (outP < 1000)

fprintf('\n      RotaryWingSim  \n');
fprintf('       u                w              q\n');
fprintf('x  %10.6f     %10.6f      %10.6f\n',xU,xW,xQ);
fprintf('z  %10.6f     %10.6f      %10.6f\n',zU,zW,zQ);
fprintf('m  %10.6f     %10.6f      %10.6f\n',mU,mW,mQ);

fprintf('\n      Measured  Data\n');
fprintf('       u                w              q\n');
fprintf('x  %10.6f     %10.6f      %10.6f\n',rfrncVal(1),rfrncVal(2),rfrncVal(3));
fprintf('z  %10.6f     %10.6f      %10.6f\n',rfrncVal(4),rfrncVal(5),rfrncVal(6));
fprintf('m  %10.6f     %10.6f      %10.6f\n\n',rfrncVal(7),rfrncVal(8),rfrncVal(9));

clcApprox(xU, xW, xQ, zU, zW, zQ, mU, mW, mQ, vErth, gloVars.mu, gloVars.alfaNf + gloVars.a1, thtaClmb, gloVars.cTs, outFile);

fprintf(outFile,'\n M _ _ M _ _ << Measured Data >> _ _ M _ _ M\n');
clcApprox(rfrncVal(1), rfrncVal(2), rfrncVal(3), rfrncVal(4), rfrncVal(5), rfrncVal(6), rfrncVal(7), rfrncVal(8), rfrncVal(9), vErth, gloVars.mu, gloVars.alfaNf + gloVars.a1, thtaClmb, gloVars.cTs, outFile);

end
end % if jOut < 2
