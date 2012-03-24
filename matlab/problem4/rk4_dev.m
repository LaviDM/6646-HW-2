function [t, y, d] = rk4_dev(f, y0, h, H, tfinal)
  % rk4_dev(f, y0, h, H, tfinal)  Uses the classical Runge-Kutta method of 
  %                               order 4 to evaluate y' = f(t, y) where y(0) 
  %                               = y0, h is the time step size, and tfinal is 
  %                               the largest t value Keeps track of the 
  %                               deviation |H(y_n) - H(y_0)|.
  
  % initialize the step number
  n = 1;
  
  % set the initial values
  y(:,n) = y0;
  t(n) = 0;
  d(n) = 0;
  
  while t(end) < (tfinal - h)
    % apply the rk4 method
    Y1 = y(:,n);
    Y2 = y(:,n) + (h/2) * f(h*n, Y1);
    Y3 = y(:,n) + (h/2) * f(h*(n+1/2), Y2);
    Y4 = y(:,n) + h * f(h*(n+1/2), Y3);
    y(:,n+1) = y(:,n) + (h/6) * (f(h*n, Y1) + 2*f(h*(n+1/2), Y2) + 2*f(h*(n+1/2), Y3) + f(h*(n+1), Y4));
    
    % calculate the deviation
    d(n+1) = norm(H(y(:,n+1)) - H(y(:,1)));
    
    % store the time step
    t(n+1) = t(n) + h;
    
    % increment the step number
    n = n + 1;
  end
end
