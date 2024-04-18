clc


mu=[0:0.01:1.0];

n= length(mu);

for i=1:length(mu)
 f(i)=(1.57 - 2*mu(i) + mu(i)*mu(i));
end
plot(mu,f)