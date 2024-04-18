function vIvHret =  lam4slv(vIvH)
% input is vIvH = (vHeli/vIhover)
%
% pSlv =  [  1/2*(-2*vIvH+2*(vIvH^2+4)^(1/2))^(1/2), 
%     -1/2*(-2*vIvH+2*(vIvH^2+4)^(1/2))^(1/2),
%   1/2*(-2*vIvH-2*(vIvH^2+4)^(1/2))^(1/2),
%  -1/2*(-2*vIvH-2*(vIvH^2+4)^(1/2))^(1/2)];


vIvHret = 1/2*(-2*vIvH+2*(vIvH^2+4)^(1/2))^(1/2);

return 

