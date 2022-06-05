  function a = coeff_InterpN_scheletro(x,y)
% Calcola la tabella delle differenze divise memorizzando solo gli
% elementi della diagonale, che rappresentano i coefficienti del polinomio
% interpolatore di Newton
n = length(x);
a=y;
for i=1:n
    for k=n-1:-1:i
        a(k+1)=(a(k+1)-a(k))/(x(k+1)-x(k+1-i)); %%
    end
end

  

