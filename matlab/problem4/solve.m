% define constants
alpha = 0.01;
beta_ = 0.6;
steps = 5000;
hsteps = [0.1, 0.01];

% Hamiltonian
H = @(y) (y(3)^2 + y(4)^2)/2 - 1/sqrt(y(1)^2 + y(2)^2) - alpha/(2*(y(1)^2 + y(2)^2)^(3/2));

% f from y' = f(t, y)
f = @(t, y) [ ...
  y(3); ...
  y(4); ...
  - y(1)*(y(1)^2 + y(2)^2)^(-3/2) - (3/2)*alpha*y(1)*(y(1)^2 + y(2)^2)^(-5/2); ...
  - y(2)*(y(1)^2 + y(2)^2)^(-3/2) - (3/2)*alpha*y(2)*(y(1)^2 + y(2)^2)^(-5/2); ...
];

% initial conditions
q1_0 = 1 - beta_;
q2_0 = 0;
p1_0 = 0;
p2_0 = sqrt((1 + beta_)/(1 - beta_));
y0 = [q1_0; q2_0; p1_0; p2_0];

p = figure;

for i=1:2
  h = hsteps(i);
  tfinal = 500;
  
  % solve via implicit midpoint
  [T, Y, D] = implicitmidpoint_dev(f, y0, h, H, tfinal);
  
  label = ['Implicit Midpoint, h = ', num2str(h,10)];
  filename = ['04_imid_', num2str(i), '.pdf'];
  plot(Y(1,:), Y(2,:));
  xlabel('y_1');
  ylabel('y_2');
  title(label);
  print(p, '-dpdf', filename);
  
  filename = ['04_imid_dev_', num2str(i), '.pdf'];
  plot(T, D);
  xlabel('time');
  ylabel('|H(q(t), p(t)) - H(q(0), p(0))|');
  title(label);
  print(p, '-dpdf', filename);
  
  % solve via rk4
  [T, Y, D] = rk4_dev(f, y0, h, H, tfinal);
  label = ['RK4, h = ', num2str(h,10)];
  filename = ['04_rk4_', num2str(i), '.pdf'];
  plot(Y(1,:), Y(2,:));
  xlabel('y_1');
  ylabel('y_2');
  title(label);
  print(p, '-dpdf', filename);
  
  filename = ['04_rk4_dev_', num2str(i), '.pdf'];
  plot(T, D);
  xlabel('time');
  ylabel('|H(q(t), p(t)) - H(q(0), p(0))|');
  title(label);
  print(p, '-dpdf', filename);
end
