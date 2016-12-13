function ALM
%%%%% Homework 4 of Numerical Optimization --Augmented Lagrange Multiplier Method

%%%%% Ch6-2-(vi) %%%%%
Ffunc=@(x1, x2) 4*x1^2+3*x2^2-5*x1*x2-8*x1;
Hfunc=@(x1, x2) x1+x2-4;
gfunc1=@(x1) -x1;
gfunc2=@(x2) -x2;


x=[1, 1];
[minX minFuncValue]=augLagMul(Ffunc, Hfunc, gfunc1, gfunc2, x, 1e-6)

end