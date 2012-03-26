function [T, Y, D] = leapfrog_dis(f, M, y0, h, H, tfinal)
  % leapfrog_dis(f, M, y0, h, H, tfinal)  
  %     Uses the leapfrog method to evaluate Mq' = p, p' = f(q) and keep track 
  %     of the Hamiltonian discrepancy where y(0) = y0, h is the time step 
  %     size, H is the Hamiltonian, and tfinal is the largest t value
  
  % initialize the step number
  n = 1;

  % set the initial values
  Y(:,n) = y0;
  T(n) = 0;
  ham0 = H(y(1,1), y(2,1));
  D(n) = H(y(1,n), y(2,n)) - ham0;
  
  % jumpstart the method
  q2(1) = Y(1,1) + (h*Y(2,1))/(2*M);
  
  while T(end) < (tfinal - h)
    % apply the symplectic method
    Y(2,n+1) = Y(2,n) + h*f(q2(n));
    q2(n+1) = (h * Y(2,n+1))/M + q2(n);
    Y(1,n+1) = (q2(n+1) + q2(n))/2;
    
    % compute the Hamiltonian discrepancy
    D(n+1) = H(Y(1,n+1), Y(2,n+1));
    
    % store the time step
    T(n+1) = T(n) + h;
    
    % increment the step number
    n = n + 1;
  end
end
