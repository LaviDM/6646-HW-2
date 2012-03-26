sigma = 10;
b = 8/3;
r = 28;

y0 = [0; 1; 0];

f = @(t, y) [sigma*(y(2) - y(1)); r*y(1) - y(2) - y(1)*y(3); y(1)*y(2) - b*y(3)];
t0 = 0;
tf = 100;
tolerance1 = 1.e-6;

options = odeset('RelTol', tolerance1);
[T, Y] = ode45(f, [t0, tf], y0, options);

p = figure;

label = 'y3 vs y1 with 1.e-6 tolerance';
filename = '05_y3_y1_e-6.pdf';
plot(Y(:,1), Y(:,3));
xlabel('y_1');
ylabel('y_3');
title(label);
print(p, '-dpdf', filename);

label = 'y2 vs t with 1.e-6 tolerance';
filename = '05_y2_t_e-6.pdf';
plot(T, Y(:,2));
xlabel('y_2');
ylabel('t');
title(label);
print(p, '-dpdf', filename);

label = '3D Lorenz with 1.e-6 tolerance';
filename = '05_3d_e-6.pdf';
plot3(Y(:,1), Y(:,2), Y(:,3));
xlabel('y_1');
ylabel('y_2');
zlabel('y_3');
title(label);
print(p, '-dpdf', filename);

y1_end = Y(end);

tolerance2 = 1.e-7;
options = odeset('RelTol', tolerance2);
[T, Y] = ode45(f, [t0, tf], y0, options);

label = 'y3 vs y1 with 1.e-7 tolerance';
filename = '05_y3_y1_e-7.pdf';
plot(Y(:,1), Y(:,3));
xlabel('y_1');
ylabel('y_3');
title(label);
print(p, '-dpdf', filename);

y2_end = Y(end);

['1.e-6 tolerance: y(100) = ', num2str(y1_end)]
['1.e-7 tolerance: y(100) = ', num2str(y2_end)]

