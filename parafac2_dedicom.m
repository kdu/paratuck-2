function [A,Bk,C,info] = parafac2_dedicom(T,R,opt)
%PARAFAC_DEDICOM PARAFAC-2 via dedicom
%  Input:
%    T - tensor IxJxK
%  Output:
%    A - IxR matrix
%    C - KxR matrix
%    Bk - JxRxK
%  so that approximately Xk = parafac2_model(A,Bk,C)
  if (~exist('opt', 'var')), opt = struct(); end;
  if (~isfield(opt, 'pow_sig')), opt.pow_sig = 1; end  

  [I,J,K] = size(T);
  %Tucker
  [Ut,St,Vt] = svd(reshape(T,I,[]));
  Ut = Ut(:,1:R);
  sVt= diag(St);

  Mc = [];%zeros(R,R,K);
  for k=1:K
    Ttemp = diag(sVt(1:R).^(-opt.pow_sig)) * Ut' *T(:,:,k);
    %Ttemp = (T(:,:,k));
    Mc(:,:,k) = Ttemp * Ttemp';
  end


  [A,SigmaEst,C, info] = pt2d_algebraic_sym(Mc,R);
  info.SigmaEst = SigmaEst;
  C = C';
  Bk = zeros(J,R,K);
 
  A = Ut * diag(sVt(1:R).^(opt.pow_sig)) * A;

  for k=1:K
    Bk(:,:,k) = ((A * diag(C(k,:))) \ T(:,:,k))';
  end
end