  function pval = pvalHornerN(a,x,xv)
% Valutazione di un polinomio espresso nella base di Newton con lo schema di Horner
% Input
%    a   vettore contenente i coefficienti del polinomio di Horner da
%    valutare
%    x:   vettore dei nodi del problema di interpolazione
%    xv: vettore dei punti in cui valutare il polinonio
%
% Output
% pval:   vettore che contiene la valutazione del polinomio di Horner in ogni
%          componente del vettore xv.
% 

  n = length(a); 
  nval=length(xv);
  pval=zeros(nval,1);
  for i=1:nval
  pval(i) = a(n);
      for k=n-1:-1:1
         pval(i)=(xv(i)-x(k))*pval(i) + a(k); %%
      end
  end
