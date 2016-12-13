function [minX minFuncValue]=augLagMul(Ffunc, Hfunc, gfunc1, gfunc2, iniGuess, tol)
LmG=ones(1,2); % two inequality constraints
LmH=1; % one equality constraints
rP=1;
gama=1.5;
rPmax=10;
x0=iniGuess;
do=true; % do-while loop
while do
    Afunc=@(x1, x2) Ffunc(x1, x2)+ ...
                    LmG(1)*max(gfunc1(x1),-LmG(1)/2/rP)+rP*max(gfunc1(x1),-LmG(1)/2/rP)+ ...
                    LmG(2)*max(gfunc2(x2),-LmG(2)/2/rP)+rP*max(gfunc2(x2),-LmG(2)/2/rP)+ ...
                    LmH*Hfunc(x1,x2)+rP*(Hfunc(x1,x2))^2;
   
    minX=powellS(Afunc, x0, 1e-6); 
    if norm(minX-x0)<tol
        do=false;
        minFuncValue=Afunc(minX(1),minX(2)); % calculate minimize finction value
    else
        x0=minX;
    end
    
    %%%% update Lamda %%%%
    LmG(1)=LmG(1)+2*rP*max(gfunc1(minX(1)),-LmG(1)/2/rP);
    LmG(2)=LmG(2)+2*rP*max(gfunc2(minX(2)),-LmG(2)/2/rP);
    LmH=LmH+2*rP*Hfunc(minX(1),minX(2));
    %%%%%%%%%%%%%%%%%%%%%%
    
    %%%% update rP %%%%
    rP=gama*rP;
    if rP > rPmax
        rP=rPmax;
    end
    %%%%%%%%%%%%%%%%%%%
end
end