function [X] = matr_sym2_decomp(Matr,I)
%matr_sym2_decomp --- convert back to partially symmetric tensor X 
% Input: 
%   Matr = nchoosek(I+1,2) x M matrix whose columns Matr(:,k)
%          are compact representations of symmetric X(:,:,k)
% Output
%   X = I x I x M tensor
  M = size(Matr,2);
  assert(size(Matr,1)==nchoosek(I+1,2));
  X = zeros(I,I,M);
  for k=1:M
    ptr = 1;
    for i=1:I
      if i > 1
        X(i,1:i-1,:) = Matr(ptr:ptr+i-2,:) / sqrt(2);
        X(1:i-1,i,:) = Matr(ptr:ptr+i-2,:) / sqrt(2);
      end
      X(i,i,:) = Matr(ptr+i-1,:);
      ptr = ptr+i;
    end
  end  
end