%---------------------------------------------------------
% This code generates Figure 3 of the following paper: 
% "Multi-Task Compressive Sensing" (Preprint, 2007)
% This example is modified from l1qc_example.m, an example 
% from l1magic.
% Coded by: Shihao Ji, ECE, Duke University
% last change: May. 15, 2007
%---------------------------------------------------------
clear all
%
base = 40; % number of initial random measurements
ns   = 100; % number of additional random measurements
load multi_results_75.mat E1_BCS E2_BCS E1_MT E2_MT;
ST_mean_75 = mean((E1_BCS+E2_BCS)/2);
ST_std_75  = std((E1_BCS+E2_BCS)/2);
MT_mean_75 = mean((E1_MT+E2_MT)/2);
MT_std_75 = std((E1_MT+E2_MT)/2);
load multi_results_50.mat E1_BCS E2_BCS E1_MT E2_MT;
ST_mean_50 = mean((E1_BCS+E2_BCS)/2);
ST_std_50  = std((E1_BCS+E2_BCS)/2);
MT_mean_50 = mean((E1_MT+E2_MT)/2);
MT_std_50 = std((E1_MT+E2_MT)/2);
load multi_results_25.mat E1_BCS E2_BCS E1_MT E2_MT;
ST_mean_25 = mean((E1_BCS+E2_BCS)/2);
ST_std_25  = std((E1_BCS+E2_BCS)/2);
MT_mean_25 = mean((E1_MT+E2_MT)/2);
MT_std_25 = std((E1_MT+E2_MT)/2);
%
ST_mean = (ST_mean_25+ST_mean_50+ST_mean_75)/3;
ST_std = (ST_std_25+ST_std_50+ST_std_75)/3;

figure; hold on;
plot((1:ns)+base,ST_mean,'k-o');
plot((1:ns)+base,MT_mean_25,'b-+');
plot((1:ns)+base,MT_mean_50,'g-s');
plot((1:ns)+base,MT_mean_75,'r-*');
xlabel('Number of Measurements'); ylabel('Reconstruction Error');
legend('ST','MT 25%','MT 50%', 'MT 75%',1); box on;
%
figure; hold on;
errorbar((1:ns)+base,ST_mean,ST_std,'k-o');
%errorbar((1:ns)+base,MT_mean_25,MT_std_25,'b-+');
%errorbar((1:ns)+base,MT_mean_50,MT_std_50,'g-s');
errorbar((1:ns)+base,MT_mean_75,MT_std_75,'r-*');
xlabel('Number of Measurements'); ylabel('Reconstruction Error');
%legend('ST','MT 25%','MT 50%', 'MT 75%',1); box on;
box on; axis([40,140,-0.2,1.5]);
legend('ST','MT 75%',1); 
