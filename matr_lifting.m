function [Tlift] = matr_lifting(Matr)
%matr_phi --- builds  tensor lifting of the matrix (tensorized Khatri-Rao)
%  Input:
%    Matr = LxK matrix
%  Output:
%    Tlift = LxLxK matrix such that Tlift(:,:,k) = Matr(:,k) * (Matr(:,k))'
 [L,K] = size(Matr);
 Tlift = zeros(L,L,K);
 for k=1:K
   Tlift(:,:,k) = Matr(:,k) * (Matr(:,k)).';
 end
end