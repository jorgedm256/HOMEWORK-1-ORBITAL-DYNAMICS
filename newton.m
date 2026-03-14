function x=newton(f,df,x,tol,maxiter)
iter=0;
for i=1:maxiter
    fvalue=f(x);
    xnew=x-f(x)/df(x);
    if abs(fvalue)<tol && (xnew-x)/(1+abs(xnew))<tol
        fprintf("The number of iterations is %.f", iter);
        return
    else
        x=xnew;
    end
    iter=iter+1;


end



if iter==maxiter
    fprintf("Convergence could not be reached");
end
end





