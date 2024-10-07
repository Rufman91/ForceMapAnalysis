function F_fit = combined_DLVO_Hertz(params, X)
    % Unpack parameters
    lambda_D = params(1); % Debye length
    epsilon = params(2); % Permittivity of medium
    H = params(3); % Hamaker constant
    sigma_1 = params(4); % Surface potential of tip
    sigma_2 = params(5); % Surface potential of substrate
    E_star = params(6); % Effective Young's modulus
    b_0 = params(7); % Offset (contact point)
    R = params(8); % Tip radius (assumed known)

    % Shift the tip height
    D = X - b_0;

    % Initialize the force array
    F_fit = zeros(size(D));

    % Calculate DLVO forces for D-b_0 > R
    DLVO_part1 = -H * R ./ (6 * D.^2);
    DLVO_part2 = (2 * pi * lambda_D * R / epsilon) * ...
                 ((sigma_1^2 + sigma_2^2) .* exp(-2 * D / lambda_D) + ...
                  2 * sigma_1 * sigma_2 .* exp(-D / lambda_D));
    
    F_DLVO = DLVO_part1 + DLVO_part2;

    % Hertzian force component where D-b_0 <= R
    F_Hertz = (4/3) * E_star * sqrt(R) .* (R - D).^(3/2);

    % Apply conditions: DLVO force for D-b_0 > R, DLVO + Hertz for D-b_0 <= R
    F_fit(D > R) = F_DLVO(D > R);
    F_fit(D <= R) = F_DLVO(D <= R) + F_Hertz(D <= R);
end
