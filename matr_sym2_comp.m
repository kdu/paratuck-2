function [Matr] = matr_sym2_comp(X)
%matr_sym2_comp ---  compact unfolding of partially symmetric tensor 
% Input: 
%   X = I x I x M tensor
% Output
%   Matr = nchoosek(I+1,2) x M matrix whose columns Matr(:,k)
%          are compact representations of (X(:,:,k)+X(:,:,k).')/2
  I = size(X,1); M = size(X,3);
  Matr = zeros(nchoosek(I+1,2), M);
  ptr = 1;
  for i=1:I
    if (i > 1)
      Matr(ptr:ptr+i-2,:) = ...
        (squeeze(X(i,1:i-1,:)) + squeeze(X(1:i-1,i,:)))/sqrt(2);
    end
    Matr(ptr+i-1,:) = X(i,i,:);
    ptr = ptr+i;
  end
end