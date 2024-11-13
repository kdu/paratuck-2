F = [1 1; 2 -1];
A = [1 1; -1 2];
B = [1 1; 1 -3];

G = [-5 -4 -3 -2 -1 0  1  2 3 4; 
      1  0  2  1 -3 2 -2 -1 0 1],
H = [-5 -4 -3 -2 -1 0 1 2 3 4; 
      1  1  1  1  1 1 1 1 1 1];
R = 2; S = 2;

T = pt2d_model(A,B,F,G,H);
[Aest,Best,Fest,Gest,Hest,info] = pt2d_algebraic_nsym(T,R,S); % Get factors
norm(T - pt2d_model(Aest,Best,Fest,Gest,Hest),'fro')^2 % Check the error

factors_dist(A,Aest)
factors_dist(B,Best)



%% Compare with ALS (ALS code needed)
addpath als

n_runs = 100;
n_success = 0;
sss_A = zeros(n_runs,1);
sss_B = zeros(n_runs,1);
figure
for i=1:n_runs
  [estT, Aest, D_A, Hest, D_B, Best, it, cost_function] = ALS_Paratuck(T,R,S);
  n_success = n_success + (cost_function(end) <=1e-10);
  semilogy(cost_function);
  sss_A(i) = factors_dist(A,Aest);
  sss_B(i) = factors_dist(B,Best);
  hold on;
end
hold off;