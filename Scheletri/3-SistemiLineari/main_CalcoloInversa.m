clc
clear all
%Calcolo dell'inversa di una matrice risolvendo n sistemi lineari aventi la
%stessa matrice; il termine noto del sistema i-esimo e' l'i-esima colonna
%della matrice identità

A1=[3 5 7; 2 3 4; 5 9 11];
A2=[1 2 3 4; 2 -4 6 8; -1 -2 -3 -1; 5 7 0 1];

%La funzione inv di matlab calcola l'inversa della matrice
invMatA1=inv(A1)

%Passiamo alla funzione sol_nsis come secondo argomento la matrice
%identità: saranno risolti n sistemi lineari aventi la
%stessa matrice; il termine noto del sistema i-esimo e' l'i-esima colonna
%della matrice identità

IA1=sol_nsis(A1,eye(3))

invMatA2=inv(A2)
IA2=sol_nsis(A2,eye(4))
