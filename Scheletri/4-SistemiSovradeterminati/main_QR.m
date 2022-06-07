clc
clear all
close all
%----------------------
x=[-3.5 -3 -2 -1.5 -0.5 0.5 1.7 2.5 3]';
y=[-3.9 -4.8 -3.3 -2.5 0.3 1.8 4 6.9 7.1]';
%x=[-3.14, -2.4, -1.57, -0.7, -0.3, 0, 0.4, 0.7, 1.57]';
% y=[0.02, -1, -0.9, -0.72, -0.2, -0.04, 0.65, 0.67, 1.1]';
%x= linspace(0, 3, 12)';
%y= exp(x).*cos(4*x) + randn(12,1);
%  x= [1.001, 1.0012, 1.0013, 1.0014, 1.0015, 1.0016]';
% y= [-1.2, -0.95, -0.9, -1.15, -1.1, -1]';
%----------------------

m=length(x);
xx=linspace(min(x),max(x),200)';


figure
hold on
plot(x,y,'k*')

n=input('grado: ');
l=length(x);

%Costruisce matrice A
A=costruisci_vander(x,n+1);

%......METODO EN................
%vettore coefficienti incogniti
[a1,condG, flag]=metodoEN(A,y);
fprintf("Indice di condizionamento A'A ---> %e \n ",condG )
% valuta approssimante
if flag~=0
disp('ATTENZIONE MATRICE EQUAZIONI NORMALI NUMERICAMENTE NON DEFINITA POSITIVA')
end
    pol1 = pvalHorner(a1,xx);


%......METODO QRLS..............
%vettore coefficienti incogniti
a2 = metodoQR(A,y); 
% valuta approssimante
pol2 = pvalHorner(a2,xx);


% Disegna il grafico dell'approssimante
 
    plot(xx,pol1,'g-',xx,pol2,'r--');
    legend('Punti','Equazioni Normali','QRLS')

%residui

%calcolo il valore del polinomio (metodo EN) nei nodi x
  pol1x=pvalHorner(a1,x);


%calcolo il valore del polinomio (metodo QRLS) nei nodi x
 pol2x=pvalHorner(a2,x);

%Calcolo la norma 2 al quadrato del residuo per ciascun metodo

res1=norm(y-pol1x);
res2=norm(y-pol2x);
fprintf("Norma 2 al quadrato del residuo metodo Equazioni Normali  -->  %e\n",res1)
fprintf("Norma 2 al quadrato del residuo metodo  --> QRLS %e",res2)
