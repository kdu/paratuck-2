%% Load data
clear all, close all;
load('test_example_4x4.mat');


T = pt2d_model(A,B,F,G,H);
I = size(T,1); J = size(T,2); K = size(T,3);
R = size(A,2); S = size(B,2);


%% Make sure we have the same linear combinations for eigenvalue computations
lc_R = [1 1; 0 1; 0 1; 0 1]; lc_S = [0 1; 1 2; 0 3; 0 1]
% lc_R =normc(rand(R,2)); lc_S = normc(rand(S,2));
opt = struct('lc_R', lc_R, 'lc_S', lc_S);



%%
%% SNR 70-30 DB 
noise_in_DB = [70 60 50 40 30];
noise_levels = sqrt(10).^(-noise_in_DB./10);

n_runs = 100;
sss_A = zeros((n_runs),length(noise_levels));
sss_B = zeros((n_runs),length(noise_levels));
mse_T = zeros((n_runs),length(noise_levels));
mse_T_hooi = zeros((n_runs),length(noise_levels));


for k=1:length(noise_levels)
  for i=1:(n_runs)
    noise = randn(I,J,K);
    noise = noise / norm(noise(:));   
    T1 = T + norm(T(:)) * noise * noise_levels(k);
  
    opt.hooi_iter = 0; % No HOOI
    [Aest1,Best1,Fest,Gest,Hest,info] = pt2d_algebraic_nsym(T1,R,S,opt); 
    sss_A(i,k) = factors_dist(A,Aest1);
    sss_B(i,k) = factors_dist(B,Best1);
    That = pt2d_model(Aest1,Best1,Fest,Gest,Hest);
    mse_T(i,k) = norm(T1-That,'fro').^2;
  
    opt.hooi_iter = 10; % HOOI with maximum 10 iterations
    [Aest1,Best1,Fest,Gest,Hest,info] = pt2d_algebraic_nsym(T1,R,S,opt); 
    sss_A_hooi(i,k) = factors_dist(A,Aest1);
    sss_B_hooi(i,k)  = factors_dist(B,Best1);
    That = pt2d_model(Aest1,Best1,Fest,Gest,Hest);
    mse_T_hooi(i,k) = norm(T1-That,'fro').^2;     
  
  end  
end
%% Plots of Noise

noise_in_DB_combined = kron(noise_in_DB,ones(1,2))

sss_A_combined(:,1:2:10) = sss_A;
sss_A_combined(:,2:2:10) = sss_A_hooi;

figure,boxplot(log10(sss_A_combined), 'Colors', ['r','b']);
xticklabels(noise_in_DB_combined);
legend(findobj(gca,'Tag','Box'),'HOOI','plain','Location','southeast');

sss_B_combined(:,1:2:10) = sss_B;
sss_B_combined(:,2:2:10) = sss_B_hooi;

figure,f2 = boxplot(log10(sss_B_combined), 'Colors', ['r','b']);
xticklabels(noise_in_DB_combined);
legend(findobj(gca,'Tag','Box'),'HOOI','plain','Location','southeast');

mse_T_combined(:,1:2:10) = mse_T/norm(T,'fro').^2;
mse_T_combined(:,2:2:10) = mse_T_hooi/norm(T,'fro').^2;
figure,f2 = boxplot(log10(mse_T_combined), 'Colors', ['r','b']);
xticklabels(noise_in_DB_combined);
legend(findobj(gca,'Tag','Box'),'HOOI','plain','Location','southeast');

%% Supporting plots of singular values
[Aest,Best,~,~,~,info] = pt2d_algebraic_nsym(T,R,S); 
figure; semilogy(info.sigPhi, 'b-o');
figure; semilogy(info.sigPA, 'b-o');
figure; semilogy(info.sigPB, 'b-o');

noise = randn(I,J,K);
noise = noise / norm(noise(:));   
T1 = T + norm(T(:)) * noise * 1e-3; %% Noise 60dB
[Aest,Best,~,~,~,info] = pt2d_algebraic_nsym(T1,R,S); 
figure; semilogy(info.sigPhi, 'b-o');
figure; semilogy(info.sigPA, 'b-o');
figure; semilogy(info.sigPB, 'b-o');
  
%%

