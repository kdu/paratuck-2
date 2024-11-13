function [F,G,H] = core_dec(C)
%core_dec --- Decomposes the Paratuck-2 core tensor
%  by an approach that uses the kernel of submatrices of Phi(T)
%  Input:
%    C - RxSxK core tensor
%  Output:
%    F - RxS matrix
%    G - RxK matrix
%    H - SxK matrix
%    so that (approximately) C_ijk = F_ij G_ik H_jk 
  [R,S,K] = size(C);

  % choose the slice that is relatively separated from 0
  F = zeros(R,S);
  F(:,1) = 1;
  F(1,:) = 1;
  
  % Compute G and H
  for i=2:R
    for j=2:S
      z1 = squeeze(C(1,1,:) .* C(i,j,:));
      z2 = squeeze(C(i,1,:) .* C(1,j,:));
      M = [z1'; z2'];
      [U,~,~] = svd(M);
      F(i,j) = -U(2,2)/U(1,2);
    end
  end

  for k=1:K
    X = C(:,:,k) ./ F;
    [U,S,V] = svd(X);
    G(:,k) = U(:,1) * sqrt(S(1,1));
    H(:,k) = V(:,1) * sqrt(S(1,1));
  end
end