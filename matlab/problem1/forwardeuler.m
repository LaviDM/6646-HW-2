function [t, y, ham] = forwardeuler(f, y0, h, H, tfinal)
  % forwardeuler(f, y0, h, H, tfinal)   Uses the forward Euler method to
  %                                     evaluate y' = f(t, y) and keep track  
  %                                     of the Hamiltonian discrepancy where 
  %                                     y(0) = y0, h is the time step size, H 
  %                                     is the Hamiltonian, and tfinal is the 
  %                                     largest t value
  
  % initialize the step number
  n = 1;

  % set the initial values
  y(:,n) = y0;
  t(n) = 0;
  ham0 = H(y(1,1), y(2,1));
  ham(n) = H(y(1,n), y(2,n)) - ham0;

  while t(end) < (tfinal - h)
    % apply the forward Euler method
    y(:,n+1) = y(:,n) + h * feval(f, t(n), y(:,n));

    % compute the Hamiltonian discrepancy
    ham(n+1) = H(y(1,n+1), y(2,n+1)) - ham0;

    % store the time step
    t(n+1) = t(n) + h;

    % increment the step number
    n = n + 1;
  end
end
