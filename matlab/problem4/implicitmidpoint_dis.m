function [T, Y, D] = implicitmidpoint_dis(f, y0, h, H, tfinal)
  % implicitmidpoint_dis(f, y0, h, H, tfinal) 
  %     Uses the implicit midpoint method to evaluate y' = f(t, y) where 
  %     y(0) = y0, h is the time step size, and tfinal is the largest t value. 
  %     Keeps track of the deviation |H(y_n) - H(y_0)|.
  
  % initialize the step number
  n = 1;
  
  % set the initial values
  Y(:,n) = y0;
  T(n) = 0;
  D(n) = 0;
  
  % create options to supress the output of fsolve
  options = optimset('Display','off');
  
  while T(end) < (tfinal - h)
    % apply the implicit midpoint method
    Y(:,n+1) = fsolve(@(x) x - Y(:,n) - h * f(h*(n+1/2), (1/2)*(x + Y(:,n))), Y(:,n), options);
    
    % store the time step
    T(n+1) = T(n) + h;
    
    % calculate the deviation
    D(n+1) = norm(H(Y(:,n+1)) - H(Y(:,1)));
    
    % increment the step number
    n = n + 1;
  end
end
