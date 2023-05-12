function peak_skewness = periodic_signal_skewness(x, y_filtered, period)
    % Find local maxima
    [maxima, maxima_indices] = findpeaks(y_filtered);

    % Calculate half-period in indices
    dx = x(2) - x(1);
    half_period_indices = round((period / 2) / dx);

    % Calculate skewness around each peak
    peak_skewness = zeros(length(maxima), 1);
    for i = 1:length(maxima)
        peak_start = max(1, maxima_indices(i) - half_period_indices);
        peak_end = min(length(y_filtered), maxima_indices(i) + half_period_indices);
        peak_signal = y_filtered(peak_start:peak_end);
        peak_skewness(i) = signal_skewness(peak_signal);
    end
end
