function [t, y, d] = implicitmidpoint_dev(f, y0, h, H, tfinal)
  % implicitmidpoint(f, y0, h, H, tfinal) Uses the implicit midpoint method to
  %                                       evaluate y' = f(t, y) where y(0) =
  %                                       y0, h is the time step size, and
  %                                       tfinal is the largest t value. Keeps
  %                                       track of the deviation |H(y_n) -
  %                                       H(y_0)|.
  
  % initialize the step number
  n = 1;
  
  % set the initial values
  y(:,n) = y0;
  t(n) = 0;
  d(n) = 0;
  
  % create options to supress the output of fsolve
  options = optimset('Display','off');
  
  while t(end) < (tfinal - h)
    % apply the implicit midpoint method
    y(:,n+1) = fsolve(@(x) x - y(:,n) - h * f(h*(n+1/2), (1/2)*(x + y(:,n))), y(:,n), options);
    
    % store the time step
    t(n+1) = t(n) + h;
    
    % calculate the deviation
    d(n+1) = norm(H(y(:,n+1)) - H(y(:,1)));
    
    % increment the step number
    n = n + 1;
  end
end
