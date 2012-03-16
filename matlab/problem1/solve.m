% define the constants
D = 90.5 * .4184e-3;
S = 1.814;
q0 = 1.41;
M = 0.9953;

% declare the initial values
q_init = 1.4155;
p_init = 1.545 * M / 48.888;

% declare the functions
% U is the scalar potential function
U = @(q) D*(1 - exp(-S*(q - q0)))^2;
% f(q) = -U_q(q)
f = @(q) -2*D*S*(1 - exp(-S*(q - q0)))*exp(-S*(q - q0));
% H is the Hamiltonian
H = @(q, p) (p^2)/(2*M) + D*(1 - exp(-S*(q - q0)))^2;
% F is the function in y' = f(t, y)
F = @(t, y) [y(2)/M; f(y(1))];

% declare the parameters
y0 = [q_init; p_init];

% possible time step sizes
hsteps = [2, 2.3684, 2.3685];

for i=1:length(hsteps)
  h = hsteps(i);

  % set the largest t value
  tfinal = h*1000;
  
  for j=1:3
    switch j
    case 1
      [t,y,ham] = forwardeuler(F, y0, h, H, tfinal);
      label = ['forward Euler, h = ', num2str(h,10)];
      filename = ['01_fe', num2str(i), '.ps'];
    case 2
      [t,y,ham] = symplecticeuler(F, y0, h, H, tfinal);
      label = ['symplectic Euler, h = ', num2str(h,10)];
      filename = ['01_se', num2str(i), '.ps'];
    case 3
      [t,y,ham] = leapfrog(f, M, y0, h, H, tfinal);
      label = ['leapfrog, h = ', num2str(h,10)];
      filename = ['01_lf', num2str(i), '.ps'];
    end
    
    p = figure;
    plot(t,ham);
    xlabel('t');
    ylabel('H(t) - H(0)');
    title(label);
    print(p, '-dps', filename);
  end
end
