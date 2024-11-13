function [T] = pt2d_model(A,B,F,G,H)
%pt2d_model  Implements the algebraic (nonsymmetric) ParaTuck-2
%  Input:
%    A,B,F,G,H = factors (IxR,JxS,RxS,RxK,SxK)
%  Output:
%    T = IxJxK tensor
%  so that
%    T(:,:,k)=A*diag(G(:,k))*F*diag(H(:,k))*B.';
  assert(size(A,2) == size(F,1));
  assert(size(B,2) == size(F,2));
  assert(size(A,2) == size(G,1));
  assert(size(B,2) == size(H,1));
  assert(size(H,2) == size(G,2));
  
  T = zeros(size(A,1), size(B,1), size(G,2));
  for k=1:size(G,2)
    T(:,:,k)=A*diag(G(:,k))*F*diag(H(:,k))*B.';
  end
end