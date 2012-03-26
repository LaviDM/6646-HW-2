alp = 0.04;
bet = 1.e+4;
gam = 3.e+7;

y0 = [1; 0; 0];
t0 = 0;

f = @(t,y) [ ...
  -alp*y(1) + bet*y(2)*y(3);  ...
  alp*y(1) - bet*y(2)*y(3) - gam*y(2)^2;  ...
  gam*y(2)^2 ...
];

tolerance = 1.e-2;
b1 = 3;
tf = b1;

options = odeset('RelTol', tolerance);
[T, Y] = ode45(f, [t0, tf], y0, options);
y1_end = Y(end,:)

b2 = 1.e+6;
tf = b2;

[T, Y] = ode15s(f, [t0, tf], y0);
y2_end = Y(end,:)