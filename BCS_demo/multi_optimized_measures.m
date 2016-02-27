%-----------------------------------------------------------------
% This code generates Figure 4 (optimized) of the following paper: 
% "Bayesian Compressive Sensing" (Preprint, 2007)
% The dataset used is similar to l1qc_example.m from l1magic package
% Coded by: Shihao Ji, ECE, Duke University
% last change: Jan. 2, 2007
%-----------------------------------------------------------------
clear all
%
total_count = 100;
N  = 512; % signal length
T  = 20;  % number of spikes
dN = 1;
base = 40; % number of initial random measurements
ns   = 80; % number of additional optimized measurements
sigma = 0.005;
%
for count = 1:total_count
    count
    rand('state', count);
    randn('state', 2*count);
    %
    % random +/- 1 signal
    x = zeros(N,1);
    q = randperm(N);
    x(q(1:T)) = sign(randn(T,1));
    % noisy observations
    A = randn(base,N);
    A = 1.01*A./repmat(sqrt(sum(A.^2,2)),[1,N]);
    e = sigma*randn(base,1);
    y = A*x + e;
    initsigma2 = std(y)^2/1e2;
    [weights,used,sigma2,errbars,basis] = BCS_fast_rvm(A,y,initsigma2,1e-8,1);

    for i = 1:ns

        K = base+i*dN;
        a = randn(dN,N);
        unused = setdiff([1:N],used);
        a(unused) = sqrt(1.01^2-1)*a(unused)/sqrt(sum(a(unused).^2)); % noise imputation
        a(used) = basis;
        % noisy observations
        e = sigma*randn(dN,1);
        t = a*x + e;
        y = [y;t];
        A = [A;a];

        initsigma2 = std(y)^2/1e2;
        [weights,used,sigma2,errbars,basis] = BCS_fast_rvm(A,y,initsigma2,1e-8,1);
        %
        xp = zeros(N,1);
        xp(used) = weights;
        err(count,i) = norm(x-xp)/norm(x);
    
    end
    
end
save optimized_results.mat err;
beep;
disp('Done!');

