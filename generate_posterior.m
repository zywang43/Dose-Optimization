function [ beta_0 , beta_1 ] = generate_posterior( s , n , u )

logitp = @(b, x) exp(b(1)+b(2).*x)./(1+exp(b(1)+b(2).*x));
prior1 = @(b1) normpdf(b1,0,100);    % prior for intercept
prior2 = @(b2) exppdf(b2,1);    % prior for slope
post = @(b) prod( binopdf( s , n ,logitp( b , u ))) ...  % likelihood
            * prior1(b(1)) * prior2(b(2));                  % priors

        
% b1 = linspace(-10, 10, 50);
% b2 = linspace(0, 10, 50);
% simpost = zeros(50,50);
% for i = 1:length(b1)
%     for j = 1:length(b2)
%         simpost(i,j) = post([b1(i), b2(j)]);
%     end;
% end;
% mesh(b2,b1,simpost)
% xlabel('Slope')
% ylabel('Intercept')
% zlabel('Posterior density')
% view(-110,30)

% cdf = @(b1,b2)integral2(@(b1,b2)post(b1,b2) , -inf , b1 , 0,b2 );
% F = rand(1);
% fun = @(b1,b2) cdf(b1,b2) - F;
% [bt] = fzero( fun , -5,  0 )

range = 5;

while 1
     b1 = ( rand(1) - 0.5 ) * range * 2;
     b2 = rand(1) * range;
     f = post( [b1 , b2] );
     r = rand(1) * 0.001;
     if r <= f
         beta_0 = b1;
         beta_1 = b2;
         break;
     end
        

end


  
