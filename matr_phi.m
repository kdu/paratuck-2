function [PhiT] = matr_phi(T)
%matr_phi --- builds a reshaped version of Phi(T) matrix from the paper
%  Input:
%    T = IxJxK tensor
%  Output:
%    PhiT = IJxIJxK matrix such that Phi(T) = reshape(PhiT, I*J*I*J,K)
  [I,J,K] = size(T);
  PhiT = matr_lifting(reshape(T, I*J, K));
end