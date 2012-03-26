% constants
alp = 0.04;
bet = 1.e+4;
gam = 3.e+7;
tolerance = 1.e-2;
b1 = 3;
b2 = 1.e+6;
t0 = 0;

% initial conditions
y0 = [1; 0; 0];

% ODE
f = @(t,y) [ ...
  -alp*y(1) + bet*y(2)*y(3);  ...
  alp*y(1) - bet*y(2)*y(3) - gam*y(2)^2;  ...
  gam*y(2)^2 ...
];

% ode options
options = odeset('RelTol', tolerance);

% solve (nonstiff)
tf = b1;
[T, Y] = ode45(f, [t0, tf], y0, options);
y1_end = Y(end,:)

% plot
p = figure;
label = 'y1 vs t (nonstiff solution)';
filename = '06_1.pdf';
plot(T, Y(:, 1));
xlabel('t');
ylabel('y1');
title(label);
print(p, '-dpdf', filename);

p = figure;
label = 'y2 vs t (nonstiff solution)';
filename = '06_2.pdf';
plot(T, Y(:, 2));
xlabel('t');
ylabel('y2');
title(label);
print(p, '-dpdf', filename);

p = figure;
label = 'y3 vs t (nonstiff solution)';
filename = '06_3.pdf';
plot(T, Y(:, 3));
xlabel('t');
ylabel('y3');
title(label);
print(p, '-dpdf', filename);

% solve (stiff)
tf = b2;
[T, Y] = ode15s(f, [t0, tf], y0);
y2_end = Y(end,:)

% plot
p = figure;
label = 'y1 vs t (semilog, stiff solution)';
filename = '06_4.pdf';
semilogx(T, Y(:, 1));
xlabel('t');
ylabel('y1');
title(label);
print(p, '-dpdf', filename);

p = figure;
label = 'y2 vs t (semilog, stiff solution)';
filename = '06_5.pdf';
semilogx(T, Y(:, 2));
xlabel('t');
ylabel('y2');
title(label);
print(p, '-dpdf', filename);

p = figure;
label = 'y2 vs t (semilog, stiff solution)';
filename = '06_6.pdf';
semilogx(T, Y(:, 3));
xlabel('t');
ylabel('y3');
title(label);
print(p, '-dpdf', filename);

% print the final y values for each solution
['nonstiff solution: y(3) = [', num2str(y1_end(1)), ', ', num2str(y1_end(2)), ', ', num2str(y1_end(3)), ']']
['stiff solution: y(1000000) = [', num2str(y2_end(1)), ', ', num2str(y2_end(2)), ', ', num2str(y2_end(3)), ']']
