function [X,FVAL] = maxllh(X,mu0,Lfunction)
ml = @(mu) (-sum(log(Lfunction(X,mu))));
[X,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN] = fminunc(ml,[mu0]);