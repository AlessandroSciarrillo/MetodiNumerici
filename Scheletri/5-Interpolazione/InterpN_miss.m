function [yv,a]= InterpN(x,y,xv)

%%%% NON presente negli scheletri (miss) %%%%

%costruzione dei coefficienti del polinomio di interpolazione di Newton
a=coeff_InterpN(x,y);
%Valutazione dell'interpolante di Newton
yv=pvalHornerN(a,x,xv);
