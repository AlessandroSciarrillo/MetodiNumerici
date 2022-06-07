clc
clear all
close all

%punti di interpolazione
L=[-55 -45 -35 -25 -15 -5 5 15 25 35 45 55 65];
T=[3.7 3.7 3.52 3.27 3.2 3.15 3.15 3.25 3.47 3.52 3.65 3.67 3.52];

% punti di valutazione per l'interpolante
xv=linspace(min(L),max(L),200);

%Interpolazione con un polinomio di grado 12 (=length(x)-1)
[yv,a]=InterpN(L,T,xv);

%-----------------------------------------------------------------------
%Valutazioni della temperatura per xx=42 con polinomio interpolante di newton
yp1=pvalHornerN(a,L,42)

%Valutazioni della temperatura per xx=-42 con polinomio interpolante di newton
yp2=pvalHornerN(a,L,-42)
%-----------------------------------------------------------------------

figure(1)
hold on
plot(xv,yv,'b-',L,T,'r*',42,yp1,'go',-42,yp2,'mo');
legend('Interpolante di Newton','Dati','stima 1','stima 2');


