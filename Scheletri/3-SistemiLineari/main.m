% test dei metodi  

clc
clear all
nmin=100;
nmax=150;
k=1;
for n=nmin:10:nmax
   A=gallery('orthog',n,1); 
   xesatta=(1:n)'; 
   b=A*xesatta;
   norm_xesatta=norm(xesatta);
   %Gauss
   %...................................
   tic
   [Lnopivot,Unopivot,flag]=LU_nopivot(A);
   Pnopivot=eye(n);
   [xnopivot,flag]=LUsolve(Lnopivot,Unopivot,Pnopivot,b);
   toc_nop(k)=toc;
   errore_relativo_nopivot(k)=norm(xesatta-xnopivot)/norm_xesatta;
   %Gauss_parziale  
   %...................................
   tic
   [Lparziale,Uparziale,Pparziale,flag]=LU_parziale(A);
   [xparziale,flag]=LUsolve(Lparziale,Uparziale,Pparziale,b);
   toc_par(k)=toc;
   errore_relativo_parziale(k)=norm(xesatta-xparziale)/norm_xesatta;
   %Matlab
   %................................... 
   tic
   xmatlab=A\b;
   toc_mat(k)=toc;
   errore_relativo_matlab(k)=norm(xesatta-xmatlab)/norm_xesatta; 
   k=k+1;
end
xn=(nmin:10:nmax)';
figure(1)
semilogy(xn,errore_relativo_nopivot,xn,errore_relativo_parziale,xn,errore_relativo_matlab)
title('Metodo di Eliminazione di Gauss')
legend('nopivot','parziale','matlab')
xlabel('Dimensione n')
ylabel('Log errore relativo')
figure(2)
curva=xn.^3;
semilogy(xn,toc_nop,xn,toc_par,xn,toc_mat,xn,curva)
title('Tempi di calcolo')
legend('nopivot','parziale','matlab','n^3')
xlabel('Dimensione n')
ylabel('Log tempo di calcolo')
