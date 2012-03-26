function [T, Y] = rk4_qr(f, y0, h, tfinal)
  % rk4_qr(f, y0, h, tfinal) 
  %     Uses the classical Runge-Kutta method of order 4 to evaluate 
  %     y' = f(t, y) where y(0) = y0, h is the time step size, and tfinal is 
  %     the largest t value. In addition, orthogonalizes y on each step.
  
  % initialize the step number
  n = 1;
  
  % set the initial values
  Y(:,n) = y0;
  T(n) = 0;
  
  while T(end) < (tfinal - h)
    % apply the rk4 method
    Y1 = Y(:,n);
    Y2 = Y(:,n) + (h/2) * f(h*n, Y1);
    Y3 = Y(:,n) + (h/2) * f(h*(n+1/2), Y2);
    Y4 = Y(:,n) + h * f(h*(n+1/2), Y3);
    Ytemp = Y(:,n) + (h/6) * (f(h*n, Y1) + 2*f(h*(n+1/2), Y2) + 2*f(h*(n+1/2), Y3) + f(h*(n+1), Y4));
    
    % perform a QR factorization to project orthogonally
    [Q, R] = qr2([Ytemp(1), Ytemp(2); Ytemp(3), Ytemp(4)]);
    Y(:,n+1) = [Q(1,1); Q(1,2); Q(2,1); Q(2,2)];
    
    % store the time step
    T(n+1) = T(n) + h;
    
    % increment the step number
    n = n + 1;
  end
end
