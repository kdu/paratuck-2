function [F,G,H] = core_dec_simple(C)
%core_dec_simple --- Decomposes the Paratuck-2 core tensor
%  by a naive algorithm that assumes that C_:,:,k \neq 0
%  Input:
%    C - RxSxK core tensor
%  Output:
%    F - RxS matrix
%    G - RxK matrix
%    H - SxK matrix
%    so that (approximately) C_ijk = F_ij G_ik H_jk 
  [R,S,K] = size(C);

  % choose the slice that is relatively separated from 0
  Creshaped = abs(reshape(C,R*S,[]));
  [~,k0] = min((max(Creshaped,[],1)./min(Creshaped,[],1)).');
  F = C(:,:,k0);
  
  % Compute G and H
  G = zeros(R,K); H = zeros(S,K);
  for k=1:K
    X = C(:,:,k) ./ F;
    [U,S,V] = svd(X);
    G(:,k) = U(:,1) * sqrt(S(1,1));
    H(:,k) = V(:,1) * sqrt(S(1,1));
  end
end