% parameters
D = 90.5 * .4184e-3;
S = 1.814;
q0 = 1.41;
M = 0.9953;

% initial values
q_init = 1.4155;
p_init = 1.545 * M / 48.888;
y0 = [q_init; p_init];

% scalar potential function
U = @(q) D*(1 - exp(-S*(q - q0)))^2;
% f(q) = -U_q(q)
f = @(q) -2*D*S*(1 - exp(-S*(q - q0)))*exp(-S*(q - q0));
% Hamiltonian
H = @(q, p) (p^2)/(2*M) + D*(1 - exp(-S*(q - q0)))^2;
% ODE
F = @(t, y) [y(2)/M; f(y(1))];

% time step sizes
hsteps = [2, 2.3684, 2.3685];

for i=1:length(hsteps)
  h = hsteps(i);

  % set the largest t value
  tfinal = h*1000;
  
  % solve using each method
  for j=1:3
    switch j
    case 1
      [T, Y, D] = forwardeuler_dis(F, y0, h, H, tfinal);
      label = ['forward Euler, h = ', num2str(h,10)];
      filename = ['01_1_', num2str(i), '.pdf'];
    case 2
      [T, Y, D] = symplecticeuler_dis(F, y0, h, H, tfinal);
      label = ['symplectic Euler, h = ', num2str(h,10)];
      filename = ['01_2_', num2str(i), '.pdf'];
    case 3
      [T, Y, D] = leapfrog_dis(f, M, y0, h, H, tfinal);
      label = ['leapfrog, h = ', num2str(h,10)];
      filename = ['01_3_', num2str(i), '.pdf'];
    end
    
    % plot and save the figures
    p = figure;
    plot(T, D);
    xlabel('t');
    ylabel('H(t) - H(0)');
    title(label);
    print(p, '-dpdf', filename);
  end
end
