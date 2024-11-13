function [T] = parafac2_model(A,Bk,C)
%parafac2_model Builds a PARAFAC-2 structured tensor
%  Input:
%    A - IxR matrix
%    C - KxR matrix
%    Bk - tensor JxRxK
%  Output:
%    T - tensor T(:,:,k) = A *diag(C(k,:)) * (Bk(:,:k))';
  K = size(C,1);
  R = size(A,2);
  J = size(Bk,1);
  I = size(A,1);
  assert(size(Bk,3) == K);

  T = zeros(I,J,K);
  for k=1:K
     T(:,:,k) = A *diag(C(k,:)) * (Bk(:,:,k))';
  end
end