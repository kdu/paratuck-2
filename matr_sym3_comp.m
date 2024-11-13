function [Matr] = matr_sym3_comp(X)
%matr_sym2_comp ---  compact unfolding of partially symmetric tensor 
% Input: 
%   X = I x I x I x M tensor
% Output
%   Matr = nchoosek(I+2,3) x M matrix whose columns Matr(:,k)
%          are compact representations of symmetrizations of 
%                  X(:,:,:,k)
  I = size(X,1); M = size(X,4);
  Matr = zeros(nchoosek(I+2,3), M);
 
  % diagonal elements:
  for i=1:I
    Matr(i,:) = X(i,i,i,:);
  end
  ptr = I+1;
  % elements (i,i,j) and (j,j,i) for (i<j)
  for i=1:I-1
    for j=i+1:I
      Matr(ptr,:) = (X(i,i,j,:) + X(i,j,i,:)+ X(j,i,i,:)) / sqrt(3);
      Matr(ptr+1,:) = (X(i,j,j,:) + X(j,i,j,:)+ X(j,j,i,:)) / sqrt(3);
      ptr = ptr+2;
    end
  end

  % elements (i,i,j) and (j,j,i) for (i<j)
  for i=1:I-2
    for j=i+1:I-1
      for k=j+1:I
        Matr(ptr,:) = (X(i,j,k,:) +X(j,i,k,:)+ X(k,i,j,:) + ...
                       X(i,k,j,:) +X(j,k,i,:)+ X(k,j,i,:)) / sqrt(6);
        ptr = ptr+1;
      end
    end  
  end
end