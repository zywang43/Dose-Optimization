function [ k_rec, percent ] = TS_mono_onepara( K , n , thre , p_real , prior )
S = zeros( 1 , K );
n_choose = zeros( 1 , K );
p_hat = zeros( 1 , K );
p = zeros( 1 , n );
I = zeros( 1 , n );
k_rec = zeros( 1, K );
u = zeros( 1 , K );
beta_1_p = 1;

u =  - 1 ./ beta_1_p .* log( 1./ prior -1 );
logitp = @(b,x) exp(b.*x)./(1+exp(b.*x));

for t = 1 : n
    t
   for i = 1 : K
       if t == 1
           p_hat( i ) = prior( i );
       else
           p_hat( i ) = logitp( beta_0 , u(i) );
       end
   end
   p_hat
   [~ ,I(t)] = min( abs( p_hat - thre ));
   p(t) = binornd( 1 ,p_real( I(t) ) );
   S( I(t) ) = S( I(t) ) + p(t);
   n_choose( I(t) ) = n_choose( I(t) ) + 1;
   [ beta_0 ] = generate_posterior_oneparalogistic( S , n_choose , u )
   
end
[ ~ , k_rec_ind ] = max( S./ n_choose );
k_rec( k_rec_ind ) = 1;
percent = n_choose ./ n;

