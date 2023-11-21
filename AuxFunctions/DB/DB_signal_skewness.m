function skewness_coeff = signal_skewness(y)
    y_mean = mean(y);
    y_std = std(y);
    
    % Compute the skewness coefficient
    skewness_coeff = mean(((y - y_mean) / y_std).^3);
end
