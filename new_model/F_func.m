function F=F_func(p,s,n,D)
if (length(p)~=length(s))||(length(p)~=length(s))
    error("The length of p, s, D should be equal but length of p = %d, s = %d, D = %d",length(p),length(s),length(D))
end

K=length(p);
F=zeros(size(p));

if n<=1
    error("n should > 1 but the input is n = %.2f",n);
    return
end

right=log(n)+log(log(n));
f=@(p,q)p.*log(p./q)+(1-p).*log((1-p)./(1-q));

for i=1:K
    if D(i)==0
        continue
    end
    if s(i)==0
        F(i)=1;
        continue
    end
    if p(i) == 0
        f=@(p,q)log(1/(1-q));
    end
    
    if p(i) == 1
        F(i)=1;
    end
    %     if s(i)*f(p(i),p(i))>=right
    %         fprintf("p = %.2f, s*f(p,p) = %.2f >= rhs = %.2f, q is an empty set.\n",p(i),s(i)*f(p(i),p(i)),right);
    %     else
    target = @(q)abs(right-s(i)*f(p(i),q));
    options = optimoptions('fmincon','OptimalityTolerance',1e-2,'Display','off');
    F(i) = fmincon(target,(1+p(i))/2,[1;-1],[1;p(i)],[],[],[],[],[],options);
    %     end
end

end