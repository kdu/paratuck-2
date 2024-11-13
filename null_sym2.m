function [Usym,sig] = null_sym2(X, nker)
%null_sym2  --- returns the symmetric nullspace of a tensor via the SVD
% Input: 
%   X = R x R x M tensor
%   nker =  number of desired vectors in the symmetric left kernel
% Output
%   Usym = R x R x nker tensor containing the symmetric nullspace
%   sig = singular values
% 1) The slices X(:,:,k) are first symmetrized
% 2) Usym(:,:,l), l=1,...,nker is an orthonormal set of symmetric matrices 
%    that are approximately orthogonal to all X(:,:,k)
  [R,S,~] = size(X);
  assert(R==S);

  % Build a compact (weighted) representation of X
  [UA,SA,~] = svd(matr_sym2_comp(X)); % diag(SA)'
  if (~exist('nker', 'var')) 
    nker = (size(UA,1) - sum(diag(SA) > 1e-8));
  end     
  sig = diag(SA);
  UA = UA(:,end-nker+1:end);
  Usym = matr_sym2_decomp(UA,R);
end