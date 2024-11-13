function [X] = matr_sym3_decomp(Matr,I)
%matr_sym2_decomp --- convert back to partially symmetric tensor X 
% Input: 
%   Matr = nchoosek(I+1,2) x M matrix whose columns Matr(:,k)
%          are compact representations of symmetric X(:,:,k)
% Output
%   X = I x I x M tensor
  M = size(Matr,2);
  assert(size(Matr,1)==nchoosek(I+2,3));
  X = zeros(I,I,I,M);
  
  % diagonal elements:
  for i=1:I
    X(i,i,i,:) = Matr(i,:);
  end
  ptr = I+1;
  % elements (i,i,j) and (j,j,i) for (i<j)
  for i=1:I-1
    for j=i+1:I
      v1 = Matr(ptr,:) / sqrt(3);
      v2 = Matr(ptr+1,:) / sqrt(3);
      X(i,i,j,:) = v1; X(i,j,i,:) = v1; X(j,i,i,:) = v1; 
      X(i,j,j,:) = v2; X(j,i,j,:) = v2; X(j,j,i,:) = v2;
      ptr = ptr+2;
    end
  end

  % elements (i,i,j) and (j,j,i) for (i<j)
  for i=1:I-2
    for j=i+1:I-1
      for k=j+1:I
        v = Matr(ptr,:) / sqrt(6);
        X(i,j,k,:) = v; X(j,i,k,:) = v; X(k,i,j,:) = v;
        X(i,k,j,:) = v; X(j,k,i,:) = v; X(k,j,i,:) = v;
        ptr = ptr+1;
      end
    end  
  end
end