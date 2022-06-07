function [x,flag] = LUsolve(L,U,P,b)
%Risolve il sistema Ax=b a partire dalla fattorizzazione
Pb=P*b
[y,flag]=Lsolve(L,Pb);
if flag==0
    [x,flag]=Usolve(U,y);
else
    return   
end
