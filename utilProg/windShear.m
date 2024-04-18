clear all
close all
clc






zM = 300;
% aM =1200;
% dM = 400;
% dltA = 30;
% dltA = deg2rad(dltA);
% alfN = 89;
% alfN = deg2rad(alfN);



alf=arcos(




z=[0:0.01:zM];

wN=30;

aLf=1/7;


for j=1:length(z)
  w(j) = wN*(z(j)/zM)^aLf;
end

plot(w,z);
grid