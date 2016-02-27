%--------------------------------------------------------------
% This code generates Figure 3 (MT 50%) of the following paper: 
% "Multi-Task Compressive Sensing" (Preprint, 2007)
% This example is modified from l1qc_example.m, an example 
% from l1magic.
% Coded by: Shihao Ji, ECE, Duke University
% last change: May. 15, 2007
%--------------------------------------------------------------
clear all
%
total_count = 100;
N = 512; % signal length
T = 20;  % number of spikes
dN = 1;
base = 40; % number of initial random measurements
ns   = 100; % number of additional random measurements
sigma = 0.005;

for count = 1:total_count
    count
    rand('state', count);
    randn('state', 2*count);
    %
    % random +/- 1 signal
    x1 = zeros(N,1);
    x2 = zeros(N,1);
    q = randperm(N);
    x1(q(1:T)) = sign(randn(T,1));
    x2(q([1:10,T+1:T+10])) = sign(randn(T,1)); % 50% similarity

    % projection matrix
    A1 = randn(base,N);
    A1 = A1./repmat(sqrt(sum(A1.^2,2)),[1,N]);
    A2 = randn(base,N);
    A2 = A2./repmat(sqrt(sum(A2.^2,2)),[1,N]);
    A{1} = A1;
    A{2} = A2;
    % noisy observations
    e1 = sigma*randn(base,1);
    e2 = sigma*randn(base,1);
    y1 = A1*x1 + e1;
    y2 = A2*x2 + e2;
    y{1} = y1;
    y{2} = y2;

    for i = 1:ns

        K = base+i*dN;
        a1 = randn(dN,N);
        a1 = a1/sqrt(sum(a1.^2));
        % noisy observations
        e1 = sigma*randn(dN,1);
        y1 = a1*x1 + e1;
        y{1} = [y{1};y1];
        A{1} = [A{1};a1];        
        %
        a2 = randn(dN,N);
        a2 = a2/sqrt(sum(a2.^2));
        % noisy observations
        e2 = sigma*randn(dN,1);
        y2 = a2*x2 + e2;
        y{2} = [y{2};y2];
        A{2} = [A{2};a2];        
        
        % solve task 1
        a = 1e2/0.1; b = 1;
        weights = mt_CS(A{1},y{1},a,b,1e-8);
        x1_BCS = weights;

        % solve task 2
        a = 1e2/0.1; b = 1;
        weights = mt_CS(A{2},y{2},a,b,1e-8);
        x2_BCS = weights;

        % solve 1&2 by MT-CS
        a = 1e2/0.1; b = 1;
        weights = mt_CS(A,y,a,b,1e-8);
        x1_MT = weights(:,1);
        x2_MT = weights(:,2);

        % reconstruction error
        E1_BCS(count,i) = norm(x1-x1_BCS)/norm(x1);
        E2_BCS(count,i) = norm(x2-x2_BCS)/norm(x2);
        E1_MT(count,i) = norm(x1-x1_MT)/norm(x1);
        E2_MT(count,i) = norm(x2-x2_MT)/norm(x2);
        
    end
    save multi_results_50.mat E1_BCS E2_BCS E1_MT E2_MT;
end
save multi_results_50.mat E1_BCS E2_BCS E1_MT E2_MT;
disp('Done!');
beep;
