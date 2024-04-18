function cT6NN = clcCTtwst(mu,lamNf,thtaN,oM,cT)
global sigma nBlade gamRot rhoAir bTl thta1...
        aLift rRot aDiskMR cGA Gmod

oMR = oM*rRot;    
    
cNNcT= rhoAir*aDiskMR*oMR^2;

rotThrst=cT*cNNcT;

[eNTwst, e1Twst, e2Twst, n1Twst, n2Twst] = clcCoeffTwst(mu, lamNf, thtaN, oMR, rotThrst);
[aN, a1, a2, b1, b2, dQ] = clcRotAnglTwst(mu, lamNf, thtaN, oMR, rotThrst);



cT6NNn = (0.5*lamNf*(bTl^2 + 0.5*mu^2) + thtaN*(bTl^3/3 + 0.5*mu^2*bTl^2 - 4/9/pi*mu^3) +...
    thta1*(0.25*bTl^4 +0.25*mu^2*bTl^2 - mu^4/32)) + mu^3*a1/8 + 0.25*mu^2*b2*bTl ;

dcT6NN = eNTwst*(0.25*(bTl^4 + mu^2*bTl^2) - mu^4/32) + n1Twst*mu*bTl^3/3 - 0.125*mu^2*e2Twst*bTl^2;

% dcT6NN = 0;
cT6NN  = (cT6NNn + dcT6NN);

