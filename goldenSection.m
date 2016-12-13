function alpha = goldenSection(func, alphaLower, alphaUpper)
%find an optimized answer with gS method.
%%%%%%
tol=1e-6;
tow=(3-sqrt(5))/2;
N=ceil(log(tol)/log(1-tow)+3);


alpha1=(1-tow)*alphaLower+tow*alphaUpper;
alpha2=tow*alphaLower+(1-tow)*alphaUpper;
funcValue1=func(alpha1);
funcValue2=func(alpha2);
for k=4:N
    if funcValue1>funcValue2
        alphaLower=alpha1;
        %funcValueLower=funcValue1;
        alpha1=alpha2;
        funcValue1=funcValue2;
        alpha2=tow*alphaLower+(1-tow)*alphaUpper;
        funcValue2=func(alpha2);
    else
        alphaUpper=alpha2;
        %funcValueUpper=funcValue2;
        alpha2=alpha1;
        funcValue2=funcValue1;
        alpha1=(1-tow)*alphaLower+tow*alphaUpper;
        funcValue1=func(alpha1);
    end
end
alpha=(alpha1+alpha2)/2;
end