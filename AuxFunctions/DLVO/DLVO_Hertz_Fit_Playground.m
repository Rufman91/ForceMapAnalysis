% Full Script for Fitting Combined DLVO and Hertzian Model with Normalization

% Step 1: Define the X (tip height) and Y (force) data
% Replace the following with your actual experimental data


% Step 2: Define the initial guess and bounds for the fitting parameters
initial_guess = [3e-9, 80*8.854e-12, 1e-20, 0.1, 0.1, 1e9, 0, 20e-9]; 
% [lambda_D, epsilon, H, sigma_1, sigma_2, E_star, b_0, R]

lb = [1e-9, 60*8.854e-12, 1e-21, 0, 0, 1e0, -10e-9, 20e-9];
ub = [10e-9, 90*8.854e-12, 1e-19, 1, 1, 1e10, 10e-9, 20e-9];

% Step 3: Normalize the X and Y data
X_min = min(X);
X_range = max(X) - X_min;
X_norm = (X - X_min) / X_range;  % Normalize X to [0, 1]

Y_min = min(Y);
Y_range = max(Y) - Y_min;
Y_norm = (Y - Y_min) / Y_range;  % Normalize Y to [0, 1]



% Step 5: Set options for fitting
options = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 10000, 'MaxIter', 10000, 'Display', 'iter');

% Step 6: Fit the model using normalized data
params_fit = lsqcurvefit(@normalized_fit, initial_guess, X_norm, Y_norm, lb, ub, options);

% Step 7: Calculate the fit in the original scale
F_fit_norm = normalized_fit(params_fit, X_norm);
F_fit = F_fit_norm * Y_range + Y_min;  % Denormalize the force back to original scale

% Step 8: Plot the results
figure;
plot(X, Y, 'bo', 'DisplayName', 'Experimental Data'); hold on;
plot(X, F_fit, 'r-', 'DisplayName', 'Fit');
xlabel('Tip Height (m)');
ylabel('Force (N)');
legend('show');
title('Fit of AFM Data with Combined DLVO and Hertzian Model');
grid on;

% Display the fitted parameters
disp('Fitted parameters:');
disp(params_fit);

% Step 4: Define the normalized fitting function
function F_norm = normalized_fit(params, X_norm)
    % Extract the physical parameters from the fit
    lambda_D = params(1);
    epsilon = params(2);
    H = params(3);
    sigma_1 = params(4);
    sigma_2 = params(5);
    E_star = params(6);
    b_0 = params(7) * X_range;  % Scale b_0 back to original units
    R = params(8);

    % Calculate the original X from normalized X
    D = X_norm * X_range + X_min - b_0;

    % Calculate DLVO and Hertzian forces in original scale
    DLVO_part1 = -H * R ./ (6 * D.^2);
    DLVO_part2 = (2 * pi * lambda_D * R / epsilon) * ...
                 ((sigma_1^2 + sigma_2^2) .* exp(-2 * D / lambda_D) + ...
                  2 * sigma_1 * sigma_2 .* exp(-D / lambda_D));
    F_DLVO = DLVO_part1 + DLVO_part2;

    F_Hertz = (4/3) * E_star * sqrt(R) .* (R - D).^(3/2);

    % Combine DLVO and Hertzian forces
    F_total = zeros(size(D));
    F_total(D > R) = F_DLVO(D > R);
    F_total(D <= R) = F_DLVO(D <= R) + F_Hertz(D <= R);

    % Normalize the output force to match the scale of Y
    F_norm = (F_total - Y_min) / Y_range;
end

