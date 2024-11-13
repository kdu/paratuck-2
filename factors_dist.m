% Distance between factors with a simple correction function
function [dist,perm] = factors_dist(A,Aest)
  C = abs((A./vecnorm(A))' * (Aest./vecnorm(Aest)));
  [perm,~] = munkres(-C); % Optimal permutation by Hungarian algorithm
  dist = sum(1-diag(C(:,perm)).^2);
end