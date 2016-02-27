%------------------------------------------------------
% This code generates Figure 2 of the following paper: 
% "Bayesian Compressive Sensing" (Preprint, 2007)
% This example is modified from l1qc_example.m, an example 
% from l1magic.
% Coded by: Shihao Ji, ECE, Duke University
% last change: Jan. 2, 2007
%------------------------------------------------------
clear all
rand('state', 1);
randn('state', 2);
%
N = 512; % signal length
T = 20;  % number of spikes
K = 100; % number of CS measurements
%
% random +/- 1 signal
x = zeros(N,1);
q = randperm(N);
x(q(1:T)) = sign(randn(T,1)); 

% projection matrix
A = randn(K,N);
A = A./repmat(sqrt(sum(A.^2,2)),[1,N]);	
% noisy observations
sigma = 0.005;
e = sigma*randn(K,1);
y = A*x + e;

% solve by BP
x0 = A'*inv(A*A')*y;
% take epsilon a little bigger than sigma*sqrt(K)
epsilon =  sigma*sqrt(K)*sqrt(1 + 2*sqrt(2)/sqrt(K));                                                                                                              
tic;
x_BP = l1qc_logbarrier(x0, A, [], y, epsilon, 1e-3);
t_BP = toc;
fprintf(1,'BP number of nonzero weights: %d\n',sum(x_BP~=0));

% solve by BCS
initsigma2 = std(y)^2/1e2;
tic;
[weights,used,sigma2,errbars] = BCS_fast_rvm(A,y,initsigma2,1e-8);
t_BCS = toc;
fprintf(1,'BCS number of nonzero weights: %d\n',length(used));
x_BCS = zeros(N,1); err = zeros(N,1);
x_BCS(used) = weights; err(used) = errbars;

% reconstruction error
E_BP = norm(x-x_BP)/norm(x);
E_BCS = norm(x-x_BCS)/norm(x);

figure
subplot(3,1,1); plot(x); axis([1 N -max(abs(x))-0.2 max(abs(x))+0.2]); title(['(a) Original Signal']);
subplot(3,1,2); plot(x_BP); axis([1 N -max(abs(x))-0.2 max(abs(x))+0.2]); title(['(b) Reconstruction with BP, K=' num2str(K)]);
subplot(3,1,3); errorbar(x_BCS,err); axis([1 N -max(abs(x))-0.2 max(abs(x))+0.2]); title(['(c) Reconstruction with BCS, K=' num2str(K)]); box on;

disp(['BP: ||I_hat-I||/||I|| = ' num2str(E_BP) ', time = ' num2str(t_BP) ' secs']);
disp(['BCS: ||I_hat-I||/||I|| = ' num2str(E_BCS) ', time = ' num2str(t_BCS) ' secs']);
