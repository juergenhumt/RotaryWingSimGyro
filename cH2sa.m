function cH = cH2sa(mu,lamNf,thtaN)
% 2*cH/(sigma*aLift)

  global dDrgN aLift
  
  
  aN = deg2rad(2.7);
  
  aN = eq1(mu,lamNf,thtaN);
  a1 = eq2(mu,lamNf,thtaN);
  b1 = eq3(mu,lamNf,thtaN);
  
  cH= mu*dDrgN/2/aLift + thtaN*a1/3 - mu*lamNf*thtaN/2 +3/4*lamNf*a1 + mu*a1^2/4- aN*b1/6 + mu*aN^2;
   
 