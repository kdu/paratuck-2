%%% Test paratuck decomposition for random decomposition of (R,S)=(4,4)
clear, clc,

% Parameters
I = 10; J = 10; 
R=4; S=4;K=300;


n_examples = 1000;
%n_runs = 100;
rel_err = zeros(n_examples,1); 
sss_A = zeros(n_examples,1);
sss_B = zeros(n_examples,1);


for j=1:n_examples  
  % random factors
  A=randn(I,R); B=randn(J,S);
  G = randn(R,K); H = randn(S,K);
  F = rand(R,S) + 0.5;

  T = pt2d_model(A,B,F,G,H);

  [Aest,Best,Fest,Gest,Hest, info] = pt2d_algebraic_nsym(T,R,S); % Get factors
      
  rel_err(j) = ...
       (norm(T - pt2d_model(Aest,Best,Fest,Gest,Hest),'fro')/ norm(T,'fro'))^2;
  sss_A(j)  = factors_dist(A,Aest);
  sss_B(j)  = factors_dist(B,Best);  
end

%%
max(rel_err)
min(rel_err)

%%

figure,boxplot (log10(rel_err));
xticklabels([]);
figure,boxplot(sss_A);
xticklabels([]);
figure,boxplot(sss_B);
xticklabels([]);