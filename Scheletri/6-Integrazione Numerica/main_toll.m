clc
clear all
format short e

syms x real 
%.......................................
 %fx=log(x); a=1; b=2; 
%fx= 3*x^3+4*x^2+3*x+1; a=-2; b=3;
  %fx= 2*x+1; a=-2; b=3;

fx=sqrt(x); a=0; b=1;
% fx=abs(x); a=-1; b=1; 
%.......................................

f=matlabFunction(fx);
I=double(int(fx,x,a,b)); 

[INt,Nt]=traptoll(f,a,b,10^(-6));
if Nt>0
fprintf('Numero intervalli con trapcomp: %d --> Valore Integrale %e \n',Nt,INt)
end
disp(' ');

[INs,Ns]=simptoll(f,a,b,10^(-6));
if Ns>0
fprintf('Numero intervalli con simpcomp: %d -> Valore Integrale %e \n',Ns,INs)
end