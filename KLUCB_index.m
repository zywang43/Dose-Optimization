function [ klucb ] = KLUCB_index( p_hat , c , n , ln)
precision = 1e-6;
lowerbound = 0;
upperbound = 1;
max_iter = 50;
d = log( ln ) + c *log( log(ln) );
value = max( p_hat , lowerbound );
u = upperbound;
count_iteration = 0;
while  count_iteration < max_iter && u - value > precision
        count_iteration = count_iteration + 1;
        m = (value + u) * 0.5;
        if n*( p_hat *log( p_hat /m ) + ( 1-p_hat )*log( (1-p_hat )/(1-m) )) > d
            u = m;
        else
            value = m;
        end
end
klucb = ( value + u ) * 0.5;

end