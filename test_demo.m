%%% Test paratuck decomposition:
%% Build input tensor
clear, clc,

% Parameters
I = 10; J = 10; 
R=4; S=4;K=300;
% random factors
A=randn(I,R); B=randn(J,S);
F=randn(R,S); G = randn(R,K); H = randn(S,K);

%% Test non-symmetric
% Build the test tensor
T = pt2d_model(A,B,F,G,H);
% Decompose with the algebraic algorithm
assert(K >= nchoosek(R+1,2) * nchoosek(S+1,2))
[Aest,Best,Fest,Gest,Hest,info] = pt2d_algebraic_nsym(T,R,S); % Get factors
norm(T - pt2d_model(Aest,Best,Fest,Gest,Hest),'fro') % Check the error
% The two matrices below should be permutation matrices
pinv(A) * Aest
pinv(B) * Best

% Plot results
figure;
semilogy(info.sigPhi);
title('Singular values of Phi')

figure;
plot(info.sigPA);
title('Singular values of P_A')
figure;
plot(info.sigPB);
title('Singular values of P_B')


%% Test symmetric ParaTuck-2 (DEDICOM)
FS = F*F'; % create a symmetric matrix
T = pt2d_model(A,A,FS,G,G);
assert(K >= nchoosek(R+3,4));
[Aest,FSest,Gest] = pt2d_algebraic_sym(T,R); % Get factorsI would 
norm(T - pt2d_model(Aest,Aest,FSest,Gest,Gest),'fro') % Check the error
pinv(A) * Aest
factors_dist(A,Aest)


