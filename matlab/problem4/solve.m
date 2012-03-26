% constants
alp = 0.01;
bet = 0.6;
steps = 5000;
hsteps = [0.1, 0.01];

% Hamiltonian
H = @(y) (y(3)^2 + y(4)^2)/2 - 1/sqrt(y(1)^2 + y(2)^2) - alp/(2*(y(1)^2 + y(2)^2)^(3/2));

% ODE
f = @(t, y) [ ...
  y(3); ...
  y(4); ...
  - y(1)*(y(1)^2 + y(2)^2)^(-3/2) - (3/2)*alp*y(1)*(y(1)^2 + y(2)^2)^(-5/2); ...
  - y(2)*(y(1)^2 + y(2)^2)^(-3/2) - (3/2)*alp*y(2)*(y(1)^2 + y(2)^2)^(-5/2); ...
];

% initial conditions
q1_0 = 1 - bet;
q2_0 = 0;
p1_0 = 0;
p2_0 = sqrt((1 + bet)/(1 - bet));
y0 = [q1_0; q2_0; p1_0; p2_0];

for i=1:2
  h = hsteps(i);
  tfinal = 500;
  
  % solve via implicit midpoint
  [T, Y, D] = implicitmidpoint_dis(f, y0, h, H, tfinal);
  
  % plot
  p = figure;
  label = ['Implicit Midpoint, h = ', num2str(h,10)];
  filename = ['04_1_', num2str(i), '.pdf'];
  plot(Y(1,:), Y(2,:));
  xlabel('y_1');
  ylabel('y_2');
  title(label);
  print(p, '-dpdf', filename);
  
  % plot the figure
  p = figure;
  filename = ['04_2_', num2str(i), '.pdf'];
  plot(T, D);
  xlabel('time');
  ylabel('|H(q(t), p(t)) - H(q(0), p(0))|');
  title(label);
  print(p, '-dpdf', filename);
  
  % solve via rk4
  [T, Y, D] = rk4_dis(f, y0, h, H, tfinal);
  
  % plot
  p = figure;
  label = ['RK4, h = ', num2str(h,10)];
  filename = ['04_3_', num2str(i), '.pdf'];
  plot(Y(1,:), Y(2,:));
  xlabel('y_1');
  ylabel('y_2');
  title(label);
  print(p, '-dpdf', filename);
  
  % plot
  p = figure;
  filename = ['04_4_', num2str(i), '.pdf'];
  plot(T, D);
  xlabel('time');
  ylabel('|H(q(t), p(t)) - H(q(0), p(0))|');
  title(label);
  print(p, '-dpdf', filename);
end
