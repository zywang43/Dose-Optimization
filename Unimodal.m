function [ k_rec , percent ] = Unimodal( K , n , thre , p_real )
p_real = 1 - abs( p_real - thre );
b = zeros( 1 , K );
l = zeros( 1 , K );
L = 1;
c = 1;
k_rec = zeros( 1 , K );
I = zeros( 1 , n );
p = zeros( 1 , n );
p_hat = zeros( 1 , K );
n_choose = zeros( 1 , K );
for t = 1 : n
    if t < K + 1
        I(t) = t;
    else
        if mod( ( l( L ) - 1 ) , 3 ) == 0
            I( t ) = L;
        else
            if L - 1 < 1
                [ ~ , ind ] = max( b( L : L + 1 ) );
                I(t) = L - 1 + ind;
            else
                if L + 1 > K
                    [ ~ , ind ] = max( b( L - 1 : L ) );
                    I(t) = L - 2 + ind;
                else
                    
                    [ ~ , ind ] = max( b( L - 1 : L + 1 ) );
                    I(t) = L - 2 + ind;
                end
            end
        end
    end
    p( t ) = binornd( 1 , p_real( I( t ) ) );
    p_hat( I(t) ) = ( p_hat( I(t) ) * n_choose( I(t) ) + p(t) ) / ( n_choose( I(t) ) + 1 );
    n_choose( I(t) ) = n_choose( I(t) ) + 1;
    for i = 1 : K
        b(i) = KLUCB_index( p_hat(i) , c, n_choose(i) , l(i) );
    end
    [~ , L ] = max( b );
    l( L ) = l( L ) + 1;
    
end
[~, k_rec_ind] = max( p_hat );
k_rec( k_rec_ind ) =1 ;
percent = n_choose ./ n;

end