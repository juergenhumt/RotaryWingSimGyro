function out=fProp(eta,vA,rProp,pAvail)
  rhoAir=1.23;

  FOM= 0.78;
  vA = vA*1.609/3.6;
  
  A=pi*rProp^2;

  out = eta - vA*FOM/(0.5*(vA + sqrt(vA^2 + 2/rhoAir/A/vA*eta*pAvail)));


