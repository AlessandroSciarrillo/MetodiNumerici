function [a,condG,p]=metodoEN(A,b)
% input A = matrice m × n, m > n
% input b = vettore termine noto di m
% output a = vettore soluzione nel senso dei minimi quadrati del sistema lineare
% sovraderterminato Ax=b, ottenuto mediante  il metodo delle Equazioni
% Normali

% Risolvi il sistema delle equazioni normali A'*A a = A'b
G=A'*A; %%
%Costruzione matrice e termine noto
condG=cond(G); %%
z=A'*b; %%

%fattorizzazione di Choleski
[L,p]=chol(G,"lower"); %%

%Soluzione del sistema triangolare inferiore
if p>0 %%
    disp('G non definita positiva') %%
    a=G\z; %%
    return %%
else %%
    b1=Lsolve(L,z); %%
    a=Usolve(L',b1); %%
end %%

