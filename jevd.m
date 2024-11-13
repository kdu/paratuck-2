function [A,V] = jevd(Ms,lc)
%JEVD Computes approximate joint eigenvalue decomposition of K  RxR matrices
%  Input:
%    Ms = RxRxK tensor
%    lc = (optional) vector of size K with coefficients of the linear comb
%  Output:
%    A = KxR matrix containing eigenvalues as rows
  [R,~,K] = size(Ms);
  if (~exist('lc', 'var')) 
    lc = rand(R,1);
  end     
  
  Mbase = reshape(reshape(Ms,R*R,[]) * lc, R, R);
  [V, ~] = eig(Mbase);
  A = zeros(K,R);
  for j=1:K
    A(j,:) = diag(V \ (Ms(:,:,j) * V));
  end
end