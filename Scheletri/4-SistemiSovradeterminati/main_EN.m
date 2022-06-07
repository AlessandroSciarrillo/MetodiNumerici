clc
clear all
close all
x=[10:.5/5:10.5]';
y=[11.0320   11.1263   11.1339   11.1339   11.1993   11.1844]';
n=4;

%punti di valutazione
xv=linspace(10,10.5,100);

disp('CASO DATI xy')
%Costruisce matrice A
A=costruisci_vander(x,n+1);

[aEN,condG,flag]=metodoEN(A,y);
fprintf("Indice di condizionamento A'A ---> %e \n ",condG )
% valuta approssimante
if flag~=0
disp('ATTENZIONE MATRICE EQUAZIONI NORMALI NUMERICAMENTE NON DEFINITA POSITIVA')
end
yvEn=pvalHorner(aEN,xv);

[aQR]=metodoQR(A,y);
yvQR=pvalHorner(aQR,xv);


%perturbazione dati
%...........................
y1=y;
x1=x;
x1(2)=x(2)+0.013;
y1(2)=y(2)-0.001;
%...........................

disp('CASO DATI PERTURBATI x1y1')
%Costruisce matrice A
Ap=costruisci_vander(x1,n+1);

[aEN_p,condG,flag]=metodoEN(Ap,y1);
fprintf("Indice di condizionamento A'A ---> %e \n ",condG )
% valuta approssimante
if flag~=0
disp('ATTENZIONE MATRICE EQUAZIONI NORMALI NUMERICAMENTE NON DEFINITA POSITIVA')
end
yvEn_p=pvalHorner(aEN_p,xv);

[aQR_p]=metodoQR(Ap,y1);
yvQR_p=pvalHorner(aQR_p,xv);

figure
hold on
axis([10 10.5 11 11.3])
plot(x,y,'bo',xv,yvEn,'r-',xv,yvQR,'b--')
plot(x1,y1,'r+',xv,yvEn_p,'m-',xv,yvQR_p,'g--')
legend('dati xy','EN xy','QR xy','dati perturbati x1y1','EN x1y1','QR x1y1')

%METODO QRLS: si vede che e' piu' stabile perche' perturbando di poco i punti, 
%la curva risulta di poco perturbata




