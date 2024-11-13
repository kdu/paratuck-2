function [U,V] = hooi_Pk(Pk,niter)
  [R,~,S,~,M] = size(Pk);
  Pk1 = reshape(permute(Pk, [2 4 1 3 5]), S, S, []); 
  Pk1comp = reshape(matr_sym2_comp(Pk1), nchoosek(S+1,2), R, R, []);
  Pk2 = reshape(permute(Pk1comp, [2 3 1 4]), R, R, []);
  Pcomp = reshape((matr_sym2_comp(Pk2)),  nchoosek(R+1,2),  nchoosek(S+1,2),[]);
  [~,~,U,V] = hooi(Pcomp,nchoosek(R,2),nchoosek(S,2),niter);
end