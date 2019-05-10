% prior1 = @(b1) ;    % prior for intercept
% prior2 = @(b2)  
post = @(b2) exppdf(b2,1);    

cdf = @(b2)integral(@(b2)post(b2) , 0,b2 );
F = rand(1);
fun = @(b2) cdf(b2) - F;
[bt] = fzero( fun , 0 )

% b1 = linspace(-10, 10, 50);
% b2 = linspace(0, 10, 50);
% simpost = zeros(50,50);
% for i = 1:length(b1)
%     for j = 1:length(b2)
%         simpost(i,j) = cdf(b1(i), b2(j));
%     end;
% end;
% mesh(b2,b1,simpost)
% xlabel('Slope')
% ylabel('Intercept')
% zlabel('Posterior density')
% view(-110,30)
% cdf(0,1)