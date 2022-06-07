clear all
close all
nin=5;
nmax=100;
nprove=length(nin:5:nmax);
tempothom=zeros(nprove,1);
tempogauss=zeros(nprove,1);

k=1;
for n=nin:5:nmax
    a=2*ones(n,1);
    b=-ones(n-1,1);
    c=-ones(n-1,1);
    d=zeros(n,1);
    d(1)=1;
    d(n)=1;
    tic
    x=thomas(a,[0;b],c,d);
    tempothom(k)=toc;
    
    
    A=diag(a)+diag(b,-1)+diag(c,1);
    tic
    [L,U,flag]=LU_nopivot(A);
    P=eye(n);
    [xg,flag]=LUsolve(L,U,P,d);
    tempogauss(k)=toc;
    k=k+1;
end
figure
semilogy(nin:5:nmax,tempothom,nin:5:nmax,tempogauss)
legend('Thomas','Gauss')