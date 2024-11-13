function [U,V,Tc,sigU,sigV] = tucker2_hosvd(T, R, S)
%tucker2_hosvd --- performs the HOSVD in first two modes
%   and  balancing in the third mode
%  Input:
%     T = I x J x K tensor
%     R,S = target paratuck ranks
%
%  Output:          
%     U = I x R matrix
%     V = J x S matrix
%  so that approximately:
%     T(:,:,k) = (U * Tc(:,:,k) * Bc')
  [I,J,K] = size(T);
  Tsc = T;
  
  % Apply the HOSVD on first two modes to the rescaled tensor
  [U,S1,~] = svd(reshape(Tsc, I, []));
  U = U(:,1:R); 
  [V,S2,~] = svd(reshape(permute(Tsc, [2 1 3]), J, []));
  V = conj(V(:,1:S)); 
  Tc = core_mult2(Tsc,U',V');
  sigU = diag(S1);
  sigV = diag(S2);
end