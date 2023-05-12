function [period, amplitude,y_filtered] = find_spatial_frequency(x, y, lower_cutoff, upper_cutoff, padding_factor,verbose)
% Ensure that the input arrays have equal length
if length(x) ~= length(y)
    error('Input arrays x and y must have equal length');
end

% Validate padding factor input
if padding_factor < 1 || floor(padding_factor) ~= padding_factor
    error('Padding factor must be an integer greater than or equal to 1');
end

% Compute the FFT and the corresponding frequency array
n = length(y);
padded_n = n * padding_factor;
y_padded = [y(:)', zeros(1, padded_n - n)]; % Add zero-padding
y_fft = fft(y_padded);
dx = (max(x) - min(x)) / (n - 1);
frequencies = linspace(0, 1/dx, padded_n);

% Filter the frequency components within the specified bandwidth
bandwidth_indices = (frequencies >= 1/upper_cutoff) & (frequencies <= 1/lower_cutoff);
y_fft_filtered = y_fft .* bandwidth_indices;

% Compute the inverse FFT to obtain the filtered signal
y_filtered = ifft(y_fft_filtered, 'symmetric');
y_filtered = y_filtered(1:n); % Remove the padding from the filtered signal

% Find the dominant spatial frequency and its amplitude
[~, max_index] = max(abs(y_fft_filtered(1:floor(padded_n/2))));
dominant_frequency = frequencies(max_index);
period = 1/dominant_frequency;
amplitude = abs(y_fft_filtered(max_index))/padded_n * 2;

if verbose
    % Plot the original and filtered signals
    figure;
    subplot(2, 2, 1);
    plot(x, y);
    title('Original signal');
    xlabel('x (m)');
    ylabel('y (m)');
    
    subplot(2, 2, 3);
    plot(x, y_filtered);
    title('Filtered signal');
    xlabel('x (m)');
    ylabel('y (m)');
    
    % Plot the spectrum before and after applying the cutoff
    subplot(2, 2, 2);
    plot(frequencies(1:floor(padded_n/2)), abs(y_fft(1:floor(padded_n/2))));
    title('Spectrum before cutoff');
    xlabel('Frequency (1/m)');
    ylabel('Amplitude');
    xlim([0 1/lower_cutoff]);
    
    subplot(2, 2, 4);
    plot(frequencies(1:floor(padded_n/2)), abs(y_fft_filtered(1:floor(padded_n/2))));
    title('Spectrum after cutoff');
    xlabel('Frequency (1/m)');
    ylabel('Amplitude');
    xlim([0 1/lower_cutoff]);
end
end
