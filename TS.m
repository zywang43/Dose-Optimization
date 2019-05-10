function [ k_rec, percent ] = TS( K , n , thre , p_real , prior )
S = zeros( 1 , K );
n_choose = zeros( 1 , K );
p_hat = zeros( 1 , K );
p = zeros( 1 , n );
I = zeros( 1 , n );
k_rec = zeros( 1, K );

for t = 1 : n
   for i = 1 : K
       if t == 1
           p_hat( i ) = prior( i );
       else
           p_hat( i ) = betarnd( S(i) + 1 , n_choose(i) - S(i) + 1);
       end
   end
   [~ ,I(t)] = min( abs( p_hat - thre ));
   p(t) = binornd( 1 ,p_real( I(t) ) );
   S( I(t) ) = S( I(t) ) + p(t);
   n_choose( I(t) ) = n_choose( I(t) ) + 1;
end
[ ~ , k_rec_ind ] = min( abs( S./ n_choose - thre ) );
k_rec( k_rec_ind ) = 1;
percent = n_choose ./ n;


end