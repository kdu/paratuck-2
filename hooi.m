function [U1,U2,U1ort,U2ort,err] = hooi(T, R1, R2, niter)
%tucker2_hosvd --- performs the HOSVD in first two modes
%   and  balancing in the third mode
%  Input:
%     T = I x J x K tensor
%     R,S = target paratuck ranks
%
%  Output:          
%     U = I x R matrix
%     V = J x S matrix
%     sc3 = vector of scalings of length K
%  so that approximately:
%     T(:,:,k) = (U * Tc(:,:,k) * Bc') * sc3(k)
  [I,J,K] = size(T);
  T1 = reshape(T, I, []);
  T2 = reshape(permute(T, [2 1 3]), J, []);
  
  % Apply the HOSVD on first two modes to the rescaled tensor
  [U1,S1,~] = svd(T1); U1ort = U1(:,R1+1:end); U1 = U1(:,1:R1); 
  [U2,S2,~] = svd(T2); U2ort = U2(:,R2+1:end); U2 = U2(:,1:R2); 

  %[U,S1,~] = svd(randn(I,I)); U = U(:,1:R); 
  %[V,S2,~] = svd(randn(J,J)); V = (V(:,1:S)); 

  err = zeros(niter+1,1);
  err(1) = norm(T - core_mult2(T,U1*U1',U2*U2'), 'fro');
  for i=1:niter
    [U1,S1,~] = svd(reshape(permute(U2' * T2, [2,1,3]), I, []));
    U1ort = U1(:,R1+1:end); U1 = U1(:,1:R1); 
    [U2,S2,~] = svd(reshape(permute(U1' * T1, [2,1,3]), J, []));
    U2ort = U2(:,R2+1:end); U2 = U2(:,1:R2); 
    err(i+1) = norm(T - core_mult2(T,U1*U1',U2*U2'), 'fro');
  end  
end