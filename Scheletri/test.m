function [x1,Xm,it]=my_newtonSys(fun,jac,x0,tolx,tolf,nmax)

matjac=jac(x0);
if det(matjac)==0
    disp('la matrice dello jacobiano calcolata nell iterato precedente non è a rango massimo')
    x1=[];
    Xm=[];
    it=[];
    return
else
    s=-matjac\fun(x0);
    it=1;
    x1=x0+s;
    fx1=fun(x1);
end

Xm=[norm(s,1)/norm(x1,1)];
while it<=nmax && norm(fx1,1)>=tolf && norm(s,1)>=tolx*norm(x1,1)
    x0=x1;
    it=it+1;
    matjac=jac(x0);
    if det(matjac)==0
        disp('la matrice dello jacobiano calcolata nell iterato precedente non è a rango massimo')
        x1=[];
        Xm=[];
        it=[];
        return
    else
        s=matjac/fun(x0);
        x1=x0-s;
        fx1=fun(x1);
        Xm=[Xm;norm(s,1)/norm(x1,1)]
    end
end

if it==nmax
    disp('Il procedimento non converge con la precisione desiderata')
end