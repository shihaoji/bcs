%---------------------------------------------------------
% This code generates Figure 2 of the following paper: 
% "Multi-Task Compressive Sensing" (Preprint, 2007)
% This example is modified from l1qc_example.m, an example 
% from l1magic.
% Coded by: Shihao Ji, ECE, Duke University
% last change: May. 15, 2007
%---------------------------------------------------------
clear all
rand('state', 1);
randn('state', 2);
%
N = 512; % signal length
T = 20;  % number of spikes
K = [90,70]; % number of CS measurements
%
% random +/- 1 signal
x1 = zeros(N,1);
x2 = zeros(N,1);
q = randperm(N);
x1(q(1:T)) = sign(randn(T,1)); 
x2(q([1:T-5,T+1:T+5])) = sign(randn(T,1)); % 75% similarity

% projection matrix
A1 = randn(K(1),N);
A1 = A1./repmat(sqrt(sum(A1.^2,2)),[1,N]);
A2 = randn(K(2),N);
A2 = A2./repmat(sqrt(sum(A2.^2,2)),[1,N]);
A{1} = A1;
A{2} = A2;
% noisy observations
sigma = 0.005;
e1 = sigma*randn(K(1),1);
e2 = sigma*randn(K(2),1);
y1 = A1*x1 + e1;
y2 = A2*x2 + e2;
y{1} = y1;
y{2} = y2;

%solve task 1 
a = 1e2/0.1; b = 1;
weights = mt_CS(A1,y1,a,b,1e-8);
fprintf(1,'Task 1 number of nonzero weights: %d\n',sum(weights~=0));
x1_BCS = weights;

% solve task 2
a = 1e2/0.1; b = 1;
weights = mt_CS(A2,y2,a,b,1e-8);
fprintf(1,'Task 2 number of nonzero weights: %d\n',sum(weights~=0));
x2_BCS = weights;

% solve 1&2 by MT-BCS
a = 1e2/0.1; b = 1;
weights = mt_CS(A,y,a,b,1e-8);
fprintf(1,'Task 1 number of nonzero weights: %d\n',sum(weights(:,1)~=0));
fprintf(1,'Task 2 number of nonzero weights: %d\n',sum(weights(:,2)~=0));
x1_MT = weights(:,1);
x2_MT = weights(:,2);

% reconstruction error
E1_BCS = norm(x1-x1_BCS)/norm(x1);
E2_BCS = norm(x2-x2_BCS)/norm(x2);
E1_MT = norm(x1-x1_MT)/norm(x1);
E2_MT = norm(x2-x2_MT)/norm(x2);

figure
subplot(3,2,1); plot(x1); axis([1 N -1.2 1.2]); title(['(a) Original Signal 1']);
subplot(3,2,2); plot(x2); axis([1 N -1.2 1.2]); title(['(b) Original Signal 2']);
subplot(3,2,3); plot(x1_BCS); axis([1 N -1.2 1.2]); title(['(c) ST Reconstruction 1, n=' num2str(K(1))]); box on;
subplot(3,2,4); plot(x2_BCS); axis([1 N -1.2 1.2]); title(['(d) ST Reconstruction 2, n=' num2str(K(2))]); box on;
subplot(3,2,5); plot(x1_MT); axis([1 N -1.2 1.2]); title(['(e) MT Reconstruction 1, n=' num2str(K(1))]); box on;
subplot(3,2,6); plot(x2_MT); axis([1 N -1.2 1.2]); title(['(f) MT Reconstruction 2, n=' num2str(K(2))]); box on;

disp(['ST1: ||I1_hat-I1||/||I1|| = ' num2str(E1_BCS)]);
disp(['MT1: ||I1_hat-I1||/||I1|| = ' num2str(E1_MT)]);
disp(['ST2: ||I2_hat-I2||/||I2|| = ' num2str(E2_BCS)]);
disp(['MT2: ||I2_hat-I2||/||I2|| = ' num2str(E2_MT)]);
