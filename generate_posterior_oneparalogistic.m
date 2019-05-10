function [ beta ] = generate_posterior_oneparalogistic( s , n , u )

logitp = @(b, x) exp(b.*x)./(1+exp(b.*x));
prior1 = @(b) normpdf(b,1,100);    % prior for intercept
post = @(b) prod( binopdf( s , n ,logitp( b , u ))) ...  % likelihood
            * prior1(b) ;                  % priors

post1=post(1)     
cdf = @(b)integral(@(b)post(b) , -inf , b  );
cdf1=cdf(1)
F = rand(1);
fun = @(b)cdf(b) - F;
[ beta ] = fzero( fun , 1 );




  
