function [t, y] = rk4qr(f, y0, h, tfinal)
  % rk4qr(f, y0, h, tfinal) Uses the classical Runge-Kutta method of order 4
  %                         to evaluate y' = f(t, y) where y(0) = y0, h is the
  %                         time step size, and tfinal is the largest t value.
  %                         In addition, orthogonalizes y on each step.
  
  % initialize the step number
  n = 1;

  % set the initial values
  y(:,n) = y0;
  t(n) = 0;

  while t(end) < (tfinal - h)
    % apply the rk4 method
    Y1 = y(:,n);
    Y2 = y(:,n) + (h/2) * feval(f, h*n, Y1);
    Y3 = y(:,n) + (h/2) * feval(f, h*(n+1/2), Y2);
    Y4 = y(:,n) + h * feval(f, h*(n+1/2), Y3);
    ytemp = y(:,n) + (h/6) * (feval(f, h*n, Y1) + 2*feval(f, h*(n+1/2), Y2) + 2*feval(f, h*(n+1/2), Y3) + feval(f, h*(n+1), Y4));
    
    % perform a QR factorization to project orthogonally
    [Q, R] = qr2([ytemp(1), ytemp(2); ytemp(3), ytemp(4)]);
    y(:,n+1) = [Q(1,1); Q(1,2); Q(2,1); Q(2,2)];
    
    % store the time step
    t(n+1) = t(n) + h;
    
    % increment the step number
    n = n + 1;
  end
end
