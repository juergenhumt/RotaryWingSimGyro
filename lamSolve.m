clc
clear all
close all

 syms tanAlf x mu cT
 solve('sqrt(x^2 + mu^2)*x/mu + cT/2 - sqrt(x^2 + mu^2)*tanAlf=0',x)

