function [T, Y, D] = symplecticeuler_dis(f, y0, h, H, tfinal)
  % symplecticeuler_dis(f, y0, h, H, tfinal)  
  %     Uses the symplectic Euler method to evaluate y' = f(t, y) and keep 
  %     track of the Hamiltonian discrepancy where y(0) = y0, h is the time 
  %     step size, H is the Hamiltonian, and tfinal the largest t value.
  
  % initialize the step number
  n = 1;

  % set the initial values
  Y(:,n) = y0;
  T(n) = 0;
  ham0 = H(Y(1,1), Y(2,1));
  D(n) = H(Y(1,n), Y(2,n)) - ham0;

  while T(end) < (tfinal - h)
    % apply the forward Euler method to the second entry
    Y(:,n+1) = Y(:,n) + h * f(T(n), Y(:,n));
    
    % store the time step
    T(n+1) = T(n) + h;
    
    % apply the backward Euler method to the first entry
    Ytemp = Y(:,n+1) + h * f(T(n+1), Y(:,n+1));
    Y(1,n+1) = Ytemp(1);
    
    % compute the Hamiltonian discrepancy
    D(n+1) = H(Y(1,n+1), Y(2,n+1));
    
    % increment the step number
    n = n + 1;
  end
end
