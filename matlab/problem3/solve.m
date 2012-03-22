% define the constants
rotations = 5;
y0 = [1; 0; 0; 1];

% declare the parameters
hsizes = [1, 0.5, 0.25];
omegas = [0.5, 1, 2];

for i = 1:length(hsizes)
  h = hsizes(i);
  
  for j = 1:length(omegas)
    omega = omegas(j);
    
    % declare the function
    F = @(t, y) [omega * y(3); omega * y(4); -omega * y(1); -omega * y(2)];
    % compute tfinal
    tfinal = (2*pi*rotations)/omega;
    
    % solve the system
    [T, Y] = rk4qr(F, y0, h, tfinal);
    
    label = ['h = ', num2str(h,10), ', omega = ', num2str(omega, 10)];
    filename = ['03_y1vy2_', num2str(i), '_', num2str(j), '.pdf'];
    
    p = figure;
    plot(Y(1,:), Y(2,:));
    xlabel('y_1');
    ylabel('y_2');
    title(label);
    print(p, '-dpdf', filename);
  end
end
