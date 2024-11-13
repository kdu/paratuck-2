function [A,F,G,info] = pt2d_algebraic_sym(T,R,opt)
%pt2d_algebraic_sym  Implements the algebraic (nonsymmetric) ParaTuck-2
%  Input:
%    T = IxIxK tensor
%  Output:
%    A,F,G = factors (IxR,RxR,RxK)
%  so that approximately
%    T = pt2d_model(A, A, F, G, G)
  if (~exist('opt', 'var')), opt = struct(); end;
  if (~isfield(opt, 'lc')), opt.lc = normc(randn(R,2)); end  

  [I,J,K] = size(T);
  if (~exist('R', 'var')), R = I; end     
  assert(I>=R); assert(I==J);

  % 2. Tucker-2 HOSVD approximation
  [U,~,Tc] = hosvd_tucker2(T, R, R);
 
  % 3. Find the symmetrized kernel of PhiT
  N = (nchoosek(R+1,3) * R)/2;
  [Pk,info.sigPhi] = matr_phi_sym22_null(Tc, N); % matrix of [p_1,...,p_N]
  
  % 4-5. Compute symmetrized kernel of Psym
  % Ptens -  RxRxRx(NR) tensor (reshaped version of R^3xNR matrix Psym)
  [Qsym,info.sigPsym] = null_sym3(reshape(Pk,R,R,R,[]),R);

  % 6. Compute the EVDs and A factor
  Ac = cpd_psym_evd(Qsym,opt.lc(:,1),opt.lc(:,2));
  % 7. Combine with Tucker-2 factor and find the core tensor 
  A = U * Ac; 
  C = core_mult2(Tc,pinv(Ac),pinv(Ac));
 
  % 8. Decompose the core tensor
  [F,G,~] = core_dec(C);
end