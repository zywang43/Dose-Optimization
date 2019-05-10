clear all 
K = 6;
n = 36;
N = 2;
thre = 0.3;
 prior = [ 0.06 0.12 0.2 0.3 0.4 0.5];
% prior = rand( 1 , K );
% p_real = [ 0.05 0.1 0.15 0.3 0.45 0.5];
p_real = [0.01 0.05 0.08 0.15 0.30 0.45];

epsilon = 0.05;
c1 = 0.8;
k_rec = zeros( 1 , K );
percent= zeros( 1 , K );
k_rec_uni = zeros( 1 , K );
k_rec_TS = zeros( 1 , K );
percent_uni = zeros( 1 , K );
percent_TS = zeros( 1 , K );
k_rec_TS_mono = zeros( 1 , K );
percent_TS_mono = zeros( 1 , K );
for i = 1 : N
    i
%     [ k_rec , percent ] = Unimodal( K , n , thre , p_real );
%     k_rec_uni = k_rec_uni + k_rec;
%     percent_uni = percent_uni + percent;
%     [ k_rec , percent ] = TS( K , n , thre , p_real , prior );
%     k_rec_TS = k_rec_TS + k_rec;
%     percent_TS = percent_TS + percent;
    [ k_rec , percent ] = TS_mono_onepara( K , n , thre, p_real , prior);
    k_rec_TS_mono = k_rec_TS_mono + k_rec
    percent_TS_mono = percent_TS_mono + percent
    
end
% 
% percent_uni = percent_uni ./ N .*100
% k_rec_uni = k_rec_uni ./ N.*100
% percent_TS = percent_TS ./ N.*100
% k_rec_TS = k_rec_TS ./ N.*100
percent_TS_mono = percent_TS_mono ./ N;
k_rec_TS_mono = k_rec_TS_mono ./ N;


    