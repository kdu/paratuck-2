function [Usym, sig] = matr_phi_sym22_null(T,nker)
%matr_phi --- builds a reshaped version of Phi(T) matrix from the paper
%  Input:
%    T = I x I x nker tensor
%  Output:
%    P = II^2 x II^2 x nker tensor containing
%     sigma_(2,2)(lker(Phi(T))) with given number of vectors in the kernel
 [I,J,K] = size(T);
 assert(I==J);

 symT = matr_sym2_comp(T);
 [nullMatr,sig] = null_sym2(matr_lifting(symT),nker); % sz1 x sz1 x nker tensor

 Usym = zeros(I*I,I*I, nker);
 for k=1:nker
   slice_sym1 = reshape(matr_sym2_decomp(nullMatr(:,:,k),I), I*I,[]);
   Usym(:,:,k) = reshape(matr_sym2_decomp(slice_sym1',I),I*I,[]);
 end
end