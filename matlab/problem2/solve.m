% constants
mu = 0.0122771171;
mu_hat = 1 - mu;
period = 17.1;
steps = [100, 1000, 10000, 20000];
tfinal = period;

% initial conditions
x1_0 = 0.994;
x2_0 = 0;
y1_0 = 0;
y2_0 = -2.00158510637908252240537862224;
y0 = [x1_0; x2_0; y1_0; y2_0];

% ODE
F = @(t, y) [ ...
  y(2); ...
  y(1) + 2*y(4) - mu_hat*(y(1) + mu)*((y(1) + mu)^2 + y(3)^2)^(-3/2) - mu*(y(1)-mu_hat)*((y(1) - mu_hat)^2 + y(3)^2)^(-3/2); ...
  y(4); ...
  y(3) - 2*y(2) - mu_hat*y(3)*((y(1) + mu)^2 + y(3)^2)^(-3/2) - mu*y(3)*((y(1) - mu_hat)^2 + y(3)^2)^(-3/2)...
];

for i = 1:length(steps)
  h = period/steps(i);
  
  % solve the system
  [T, Y] = rk4(F, y0, h, tfinal);
  
  label = ['RK4, h = ', num2str(h,10)];
  filename = ['02_', num2str(i), '.pdf'];
  
  % plot the figure
  p = figure;
  plot(Y(1,:), Y(3,:));
  xlabel('u_1');
  ylabel('u_2');
  title(label);
  print(p, '-dpdf', filename);
end
