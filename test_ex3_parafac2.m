%%%% Example of PARAFAC-2 decomposition via DEDICOM
clc, clear all, close all;
% Parameters, noiseless tensor generation. Saved in mat files
%R = 3; I = 30; J = 40; K = 20;
%A = (rand(I,R) + 1);
%B = (rand(J,R)); C = (rand(K,R));%
load('test_parafac2.mat'); 

%%
[I,R] = size(A);
J = size(B,1);
K = size(C,1);

Sigma = B' * B;
%  Create B factors from circular shift
B_factors = zeros(J,R,K);
for k=1:K, B_factors(:,:,k) = circshift(B,k-1); end
% Generate the true tensor
[T0] = parafac2_model(A,B_factors,C);

tic
[Aest0,Bkest0,Cest0,info0] = parafac2_dedicom(T0,R); % Get factors
toc
figure, semilogy(info0.sigPhi,'b-o')
factors_dist(Aest0,A)

% noiseless covariance
M0 = pt2d_model(Aest0,Aest0,info0.SigmaEst,Cest0',Cest0'); 


%% Noisy test
lc_evd = [1 1; 0 1; 0 1]; % fix linear combinations for joint EVD

noise_in_DB = [50 40 30 20];
noise_levels = sqrt(10).^(-noise_in_DB./10);

n_runs = 100;
sss_A = zeros((n_runs),length(noise_levels));
mse_T = zeros((n_runs),length(noise_levels));
mse_T0 = zeros((n_runs),length(noise_levels));
mse_M0 = zeros((n_runs),length(noise_levels));


for k=1:length(noise_levels)
  for i=1:(n_runs)
     noise = randn(I,J,K); noise = noise / norm(noise(:));
     noise = noise * noise_levels(k) * norm(T0(:));
     T = T0 + noise;
     [Aest,Bkest,Cest, info] = parafac2_dedicom(T,R,struct('lc','lc_evd')); % Get factors
     
     sss_A(i,k) = factors_dist(Aest,A);
     mse_M0(i,k) = norm(M0- pt2d_model(Aest,Aest,info.SigmaEst,Cest',Cest'),'fro').^2;
     mse_T(i,k) = norm(T-parafac2_model(Aest,Bkest,Cest),'fro').^2;
     mse_T0(i,k) = norm(T0-parafac2_model(Aest,Bkest,Cest),'fro').^2;     
  end
end

%% Plots
figure;
boxplot(log10(sss_A));
xticklabels(noise_in_DB);

%%
figure;
boxplot(log10(mse_M0/ norm(M0,'fro').^2));
xticklabels(noise_in_DB);

%%
rmse_M0 = sqrt(sum(mse_M0,1)/n_runs) / norm(M0,'fro');
figure
semilogy(rmse_T0','b-s'); hold off;
xticks(1:4)
xticklabels(noise_in_DB);

%%
rmse_T = sqrt(sum(mse_T,1)/n_runs) / norm(T0,'fro');
rmse_T0 = sqrt(sum(mse_T0,1)/n_runs) / norm(T0,'fro');
figure
semilogy(rmse_T','r-o'); hold on;
semilogy(rmse_T0','b-s'); hold off;
xticks(1:4)
xticklabels(noise_in_DB);
legend('STD', 'RMSE', 'Location','northwest')




