function [T] = core_mult2(C,A,B)
%core_mult  --- Multiply core tensor by two matrices
  assert(size(C,1)==size(A,2));
  assert(size(C,2)==size(B,2));
  T = zeros(size(A,1),size(B,1),size(C,3));
  for k=1:size(C,3)
    T(:,:,k) = A * C(:,:,k) * B.';
  end
end