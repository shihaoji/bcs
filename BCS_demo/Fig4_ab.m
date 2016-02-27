%------------------------------------------------------
% This code generates Figures 4(a,b) of the following paper: 
% "Bayesian Compressive Sensing" (Preprint, 2007)
% This example is modified from l1qc_example.m, an example 
% from l1magic. 
% Coded by: Shihao Ji, ECE, Duke University
% last change: Jan. 2, 2007
%------------------------------------------------------
clear all
%
load random_results.mat
rand_mean = mean(err);
rand_std = std(err);
load optimized_results.mat
opt_mean = mean(err);
opt_std = std(err);
load approx_results.mat
app_mean = mean(err);
app_std = std(err);
%
base = 40;
ns = 80;
dN = 1;
% plot the mean
figure
plot(base+(1:ns)*dN,rand_mean,'b-o');
hold on;
plot(base+(1:ns)*dN,opt_mean,'r-*');
hold on;
plot(base+(1:ns)*dN,app_mean,'k-s');
xlabel('Number of Measurements'); ylabel('Reconstruction Error');
box on;
legend('Random','Optimized','Approx.',1);
% plot the variance
figure
errorbar(base+(1:ns)*dN,rand_mean,rand_std,'b-o');
hold on;
errorbar(base+(1:ns)*dN,opt_mean,opt_std,'r-*');
% hold on;
% errorbar(base+(1:ns)*dN,app_mean,app_std,'k-s');
xlabel('Number of Measurements'); ylabel('Reconstruction Error');
box on; axis([39,121,-0.3,1.6]);
legend('Random','Optimized',1);
