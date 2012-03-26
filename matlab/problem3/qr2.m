function [Q, R] = qr2(A)
  % qr2(A)  
  %     Computes a QR factorization of A using the Gram-Schmidt method.
  
  [m, n] = size(A);
  
  % compute QR using Gram-Schmidt
  for j = 1:n
    % jth column of A
    v = A(:,j);
    
    for i=1:j-1
      % compute jth column of R
      R(i,j) = Q(:,i)'*A(:,j);
      % subtract off parallel components
      v = v - R(i,j)*Q(:,i);
    end
    
    % compute length of v
    R(j,j) = norm(v);
    % normalize jth column vector
    Q(:,j) = v/R(j,j);
  end
end