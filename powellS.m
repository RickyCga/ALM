function X=powellS(func, iniGuess, tol)
%find an optimized answer with powell's method.
count=1; % for counting iteraction times
alpha=[-1, 1]; % search area in goldenSection from 1 to -1
S=eye(2); % two direction for searching, and initial direction is [0, 1] and [1, 0]
X=iniGuess; %load start point
Y=X; %record start point
n=size(S, 2); % two direction, so do twice for updating direction.
do=true; % create do-while loop in matlab
while do
    for q=1:n
        func_gS=@(alphaValue) func(X(1)+alphaValue*S(q,1), X(2)+alphaValue*S(q,2));
        alphaOpti=goldenSection(func_gS, alpha(1), alpha(2));
        X=X+alphaOpti*S(q,:);
        if (q==n) % after update the point, then update the direction and check whether the function converged.
            count=count+1;
            S(q+1,:)=(X-Y)/norm(X-Y); % calculate new direction
            S(1,:)=[];
            Y=X;
            alphaOpti=goldenSection(func_gS, alpha(1), alpha(2));
            X=X+alphaOpti*S(q,:);
            if norm(Y-X)<tol % check whther it converged
                do=false; % brake from while
            end
        end
        if count==100
            S=eye(2);
            count=1;
            break;
        end
    end
end
end
