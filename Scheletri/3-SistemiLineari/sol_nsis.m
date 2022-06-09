function X=sol_nsis(A,B)
% Funzione per risolvere l'equazione matriciale AX=B

%Input:
% A, B matrici

%Output:
% X matrice soluzione del sistema AX=B
n=size(B,1);
X=zeros(size(A));
Y=zeros(size(A));

[L,U,P,flag]=LU_parziale(A);

%[L,U,flag]=LU_nopivot(A);
%P=eye(n);

if flag==0
      for i=1:n
        Y(:,i)=Lsolve(L,P*B(:,i));  
        X(:,i)=Usolve(U,Y(:,i));     
      end
      
else
    disp('non posso proseguire');
    X=[];
    return
end


%%%%%%% OPPURE %%%%%%%%%
function A1=inversa(A)

I=eye(size(A));
[L,U,P,flag]=LU_nopivot(A);
A1=zeros(size(A));

for i=1:n
    b=I(:,i);
    [x,flag]=LU_SOLVE(L,U,P,b);  %Se si una gaus con pivoting b ->  P*b ;
    A1(:,i)=x;
end
