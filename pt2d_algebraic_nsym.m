function [A,B,F,G,H,info] = pt2d_algebraic_nsym(T,R,S,opt)
%pt2d_algebraic_nsym  Implements the algebraic (nonsymmetric) ParaTuck-2
%  Input:
%    T = IxJxK tensor
%  Output:
%    A,B,F,G,H = factors (IxR,JxS,RxS,RxK,SxK)
%  so that approximately
%    T = pt2d_model(A, B, F, G, H)
  [I,J,K] = size(T);
  if (~exist('R', 'var')), R = I; end     
  if (~exist('S', 'var')), S = J; end
  if (~exist('opt', 'var')), opt = struct(); end;
  if (~isfield(opt, 'lc_R')), opt.lc_R = normc(randn(R,2)); end  
  if (~isfield(opt, 'lc_S')), opt.lc_S = normc(randn(S,2)); end   
  if (~isfield(opt, 'hooi_iter')), opt.hooi_iter = 0; end   

  assert(I >= R); assert(J >= S);

  info = struct();

  % 2. Tucker-2 HOSVD approximation
  [U,V,Tc,info.sigU,info.sigV] = hosvd_tucker2(T, R, S);
 
  % 3. Find the symmetrized kernel of PhiT
  PhiTens = matr_phi(Tc); % Build RSxRSxK tensor (reshaping of Phi(T))
  [Pk,info.sigPhi] = null_sym2(PhiTens, nchoosek(R,2)*nchoosek(S,2));
  Pk = reshape(Pk, R,S,R,S,[]); % matrix of [p_1,...,p_M]
  
  % 4-5. Compute symmetrized kernel of PA and PB
  if opt.hooi_iter == 0
    % PAtens -  RxRx(MS^2) tensor (reshaped version of R^2xMS^2 matrix PA)
    % PBtens -  SxSx(MR^2) tensor (reshaped version of S^2xMR^2 matrix PB)
    PAtens = reshape(permute(Pk, [1 3 2 4 5]), R,R, []);
    PBtens = reshape(permute(Pk, [2 4 1 3 5]), S,S, []);
    [QA,info.sigPA] = null_sym2(PAtens, R);
    [QB,info.sigPB] = null_sym2(PBtens, R);
  else
   [UA,UB] = hooi_Pk(Pk,opt.hooi_iter);
    QA = matr_sym2_decomp(UA,R);
    QB = matr_sym2_decomp(UB,R);
  end



  % 4-5. Compute symmetrized kernel of PB
 
  % 6. Compute the EVDs and the core tensor
  Ac = cpd_psym_evd(QA,opt.lc_R(:,1),opt.lc_R(:,2));
  Bc = cpd_psym_evd(QB,opt.lc_S(:,1),opt.lc_S(:,2));

  % 7. Combine with Tucker-2 factor
  A = U * Ac; B = V * Bc; 

  % 8. Find the core tensor
  C = core_mult2(Tc,pinv(Ac),pinv(Bc));
 
  % 9. Decompose the core tensor
  [F,G,H] = core_dec(C);
 end