%% Initialization
K=6;        %the number of dose levels
n=36;       %the number of clinical trials
a0=1/2;     % initialize a
a_max=10000;
delta = 1/n;
prior=[ 0.06 0.12 0.2 0.3 0.4 0.5]; % the prior
p_real = [0.1 0.25 0.4 0.5 0.65 0.75];
q_real = [0.1 0.25 0.4 0.5 0.65 0.75];
theta=0.3; %MTD
%calculate dose level
d_func = @(x)1/2.*log((1+x)./(1-x)); 
d=d_func((prior.^(1/a0).*2.-1));
toxicity=@(x,a)((tanh(x)+1)/2).^a;
t=1; %trial
N=zeros(n,K); %number of selection
p_hat=zeros(n,K); %estimated toxicity
q_hat=zeros(n,K); %estimated efficacy 
a_hat=zeros(n+1,1); %estimated overall a
ak_hat=ones(n,K)*a0; % estimated individual a
w=zeros(n,K); % weight
%% Experiment
while t<=n
    t
    if t==1
        a_hat(1)=a0;
    else 
        a_hat(t)=sum(w(t-1,:).*ak_hat(t-1,:));
    end
    a=a_hat(t)
    if t==1
        [~,I]=max(prior.*(prior<=theta));%when t=1, choose the dose according to prior
        X=binornd(1,q_real(I));
        Y=binornd(1,p_real(I));
        q_hat(t,I)=X;
        p_hat(t,I)=Y;
        N(t,I)=1;
    else
        alpha=alpha_func(d,K,delta,t)% calculate alpha
        tox=toxicity(d,a_hat(t)-alpha)
        D=toxicity(d,a_hat(t)-alpha)<=theta%available set
        f=F_func(q_hat(t,:),N(t,:),t,D);
        [~,I]=max(f);
        I=I(1);
        X=binornd(1,q_real(I));
        Y=binornd(1,p_real(I));
        q_hat(t,:)=q_hat(t-1,:);
        q_hat(t,I)=(q_hat(t-1,I)*N(t-1,I)+X)/(N(t-1,I)+1);
        p_hat(t,:)=p_hat(t-1,:);
        p_hat(t,I)=(p_hat(t-1,I)*N(t-1,I)+Y)/(N(t-1,I)+1);%update q_hat, p_hat
        N(t,:)=N(t-1,:);
        N(t,I)=N(t-1,I)+1;
        ak_hat(t,:)=ak_hat(t-1,:);
    end
        ak_hat(t,I)=log(p_hat(t,I))/log((tanh(d(I))+1)/2);
        if ak_hat(t,I)>a_max
            ak_hat(t,I)=a_max;
        end
        w(t,:)=N(t,:)/t;
        t=t+1;
end
%% Output
a_hat(t)=sum(w(t-1,:).*ak_hat(t-1,:));
p_out=toxicity(d,a_hat(t));
[~,I] = max(p_out.*(p_out<=theta));
MTD=zeros(1,K);
MTD(I)=1
bar(N(n,:))
title("Recommended")
fprintf('toxic doses: %d\n', sum(N(n,:).*(p_real>theta)))




