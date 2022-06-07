clc
clear all
close all
%Definizione simbolica della funzione f, al fine di valutare l'integrale
%definito esatto mediante la funzione specifica di matlab, int.
syms x real
%...............................
%fx=x.^10; a=0; b=1;
% fx=asin(x); a=0; b=1; 
 fx=log(1+x); a=0; b=1;
%...............................
%Calcolo dell'integrale esatto con la funzione int di matlab che lavora su
%variabili simboliche e cast a double del risultato in output
Iesatto=double(int(fx,x,a,b));


%Conversione della funzione simbolica in matlab function
f=matlabFunction(fx);

i=0;
for N=[1 2 4 8 16 32 64 128 256]
    i=i+1;
    [INt1(i)]=TrapComp(f,a,b,N);
    
    [INs1(i)]=SimpComp(f,a,b,N);
    %Calcolo dell'errore relativo rispetto alla soluzione esatta (Trapezi)
    Et(i)=abs(INt1(i)-Iesatto)/abs(Iesatto);
    %Calcolo dell'errore relativo rispetto alla soluzione esatta (Simpson)
    Es(i)=abs(INs1(i)-Iesatto)/abs(Iesatto);
end
 
Et.'

Es.'
N=[1 2 4 8 16 32 64 128 256];
%...............................
figure(1)
semilogy(N,Et,'bo-',N,Es,'ro-')
legend('Trapezi','Simpson')
 
