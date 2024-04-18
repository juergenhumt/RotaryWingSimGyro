function dQmr = oMitr(  oM, vBd, omBd,thtaN, A1S, B1S, TMR)
[lamNfOut, vCt, omCt, anglC] = MainRotorCntrlFunc(vBd, omBd, oM, thtaN, A1S, B1S, TMR);
dQmr = anglC(14);
return 
