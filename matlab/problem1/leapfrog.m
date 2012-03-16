function [t, y, ham] = leapfrog(f, M, y0, h, H, tfinal)
  % leapfrog(f, M, y0, h, H, tfinal)  Uses the leapfrog method to evaluate 
  %                                   Mq' = p, p' = f(q) and keep track of the 
  %                                   Hamiltonian discrepancy where y(0) = y0, 
  %                                   h is the time step size, H is the 
  %                                   Hamiltonian, and tfinal is the largest t 
  %                                   value
  
  % initialize the step number
  n = 1;

  % set the initial values
  y(:,n) = y0;
  t(n) = 0;
  ham0 = H(y(1,1), y(2,1));
  ham(n) = H(y(1,n), y(2,n)) - ham0;
  
  % jumpstart the method
  q2(1) = y(1,1) + (h*y(2,1))/(2*M);
  
  while t(end) < (tfinal - h)
    % apply the symplectic method
    y(2,n+1) = y(2,n) + h*feval(f, q2(n));
    q2(n+1) = (h * y(2,n+1))/M + q2(n);
    y(1,n+1) = (q2(n+1) + q2(n))/2;
    
    % compute the Hamiltonian discrepancy
    ham(n+1) = H(y(1,n+1), y(2,n+1));
    
    % store the time step
    t(n+1) = t(n) + h;
    
    % increment the step number
    n = n + 1;
  end
end
