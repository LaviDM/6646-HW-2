function [T, Y, D] = forwardeuler_dis(f, y0, h, H, tfinal)
  % forwardeuler_dis(f, y0, h, H, tfinal)
  %     Uses the forward Euler method to evaluate y' = f(t, y) and keep track
  %     of the Hamiltonian discrepancy where y(0) = y0, h is the time step 
  %     size, H is the Hamiltonian, and tfinal is the largest t value.
  
  % initialize the step number
  n = 1;
  
  % set the initial values
  Y(:,n) = y0;
  T(n) = 0;
  ham0 = H(Y(1,1), Y(2,1));
  D(n) = H(Y(1,n), Y(2,n)) - ham0;
  
  while T(end) < (tfinal - h)
    % apply the forward Euler method
    Y(:,n+1) = Y(:,n) + h*f(T(n), Y(:,n));
    
    % compute the Hamiltonian discrepancy
    D(n+1) = H(Y(1,n+1), Y(2,n+1)) - ham0;
    
    % store the time step
    T(n+1) = T(n) + h;
    
    % increment the step number
    n = n + 1;
  end
end
