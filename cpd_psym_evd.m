function [A] = cpd_psym_evd(Mtens, lc1, lc2)
%cpd_psym_evd  --- returns the partiallysymmetric CPD of a tensor via the EVD
% Input: 
%   X = R x R x ... x R tensor of order d \ge 3
%   lc1 = linear combination (the base matrix)
%   lc2 = linear combination (for computing eigenvectors)
% Output
%   A = the 1-(d-1) mode factors, so that: Mtens = [[A,...,A,C]]
  R = size(Mtens,1);
  d = ndims(Mtens);

  %% Compute shift matrices
  if (~exist('lc1', 'var')), lc1 = rand(R,1); end  
  if (~exist('lc2', 'var')), lc2 = rand(R,1); end   

  MtensPerm = permute(reshape(Mtens, R,[],R), [2 3 1]);
  Mshift = zeros(R,R,R);
  M0 = reshape(reshape(MtensPerm,[], R) * lc1, [], R);    % Base matrix
  for j=1:R
    Mshift(:,:,j) = M0 \ MtensPerm(:,:,j);   
  end

  [A, ~] = jevd(Mshift,lc2);  % Joint EVD
end

