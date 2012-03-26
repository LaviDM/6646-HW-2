% constants
sigma = 10;
b = 8/3;
r = 28;
t0 = 0;
tf = 100;
tolerance1 = 1.e-6;
tolerance2 = 1.e-7;

% initial conditions
y0 = [0; 1; 0];

% ode
f = @(t, y) [sigma*(y(2) - y(1)); r*y(1) - y(2) - y(1)*y(3); y(1)*y(2) - b*y(3)];

% ode options
options = odeset('AbsTol', tolerance1);

% solve w/ tolerance1
[T, Y] = ode45(f, [t0, tf], y0, options);

% plot y3 vs y1
p = figure;
label = 'y3 vs y1 with 1.e-6 tolerance';
filename = '05_1.pdf';
plot(Y(:,1), Y(:,3));
xlabel('y_1');
ylabel('y_3');
title(label);
print(p, '-dpdf', filename);

% plot y2 vs t
p = figure;
label = 'y2 vs t with 1.e-6 tolerance';
filename = '05_2.pdf';
plot(T, Y(:,2));
xlabel('y_2');
ylabel('t');
title(label);
print(p, '-dpdf', filename);

% plot 3d
p = figure;
label = '3D Lorenz with 1.e-6 tolerance';
filename = '05_3.pdf';
plot3(Y(:,1), Y(:,2), Y(:,3));
xlabel('y_1');
ylabel('y_2');
zlabel('y_3');
title(label);
print(p, '-dpdf', filename);

y1_end = Y(end,:);

% ode options
options = odeset('AbsTol', tolerance2);

% solve w/ tolerance 2
[T, Y] = ode45(f, [t0, tf], y0, options);

% plot y3 vs y1
p = figure;
label = 'y3 vs y1 with 1.e-7 tolerance';
filename = '05_4.pdf';
plot(Y(:,1), Y(:,3));
xlabel('y_1');
ylabel('y_3');
title(label);
print(p, '-dpdf', filename);

y2_end = Y(end,:);

% print difference in final values
['1.e-6 tolerance: y(100) = [', num2str(y1_end(1)), ', ', num2str(y1_end(2)), ', ', num2str(y1_end(3)), ']']
['1.e-7 tolerance: y(100) = [', num2str(y2_end(1)), ', ', num2str(y2_end(2)), ', ', num2str(y2_end(3)), ']']

