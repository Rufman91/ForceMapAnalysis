% Sample code for generating a skewed periodic signal
clear all;

ST_Skew = linspace(0,1,1000);
figure
for j=[0]
for i=1:length(ST_Skew)
% Parameters
N = 1000; % number of data points
x = linspace(0, 10, N); % x-axis values
noise_amplitude = 0.1; % amplitude of noise

% Generate the sine wave
sine_freq = 2; % frequency of the sine wave
sine_amplitude = 1; % amplitude of the sine wave
sine_wave = sine_amplitude * sin(2 * pi * sine_freq * x);

% Generate the sawtooth wave
sawtooth_freq = 2*2*pi; % frequency of the sawtooth wave
sawtooth_amplitude = 0.7; % amplitude of the sawtooth wave
sawtooth_wave = 1 + sawtooth_amplitude * sawtooth(x * sawtooth_freq, ST_Skew(i));

slant_signal = x.*(j/max(x));

% figure
% plot(x,sine_wave,x,sawtooth_wave)

% Combine the sine and sawtooth waves
skewed_signal = sine_wave.*sawtooth_wave + slant_signal;

% Add noise to the signal
noisy_signal = skewed_signal + noise_amplitude * randn(1, N);

% Analyze the signal using the functions provided earlier
lower_cutoff = .1;
upper_cutoff = 100;
padding_factor = 4;
noisy_signal_final = noisy_signal;
noisy_signal_final = detrend(noisy_signal,3);
[period, amplitude, y_filtered] = find_spatial_frequency(x, noisy_signal_final, lower_cutoff, upper_cutoff, padding_factor,false);
peak_skewness_plus = periodic_signal_skewness(x, y_filtered, period);
peak_skewness_minus = -periodic_signal_skewness(x, -y_filtered, period);

Skew(i) = median([peak_skewness_plus ; peak_skewness_minus]);

end
subplot(2,1,1)
plot(x,noisy_signal)
hold on
legend
subplot(2,1,2)
plot(Skew')
legend
hold on
plot(range(Skew)/2*sign(Skew'))
end
% % Plot the original and filtered signals
% figure;
% subplot(2, 1, 1);
% plot(x, noisy_signal);
% title('Noisy skewed periodic signal');
% xlabel('x');
% ylabel('y');
% 
% subplot(2, 1, 2);
% plot(x, y_filtered);
% title('Filtered signal');
% xlabel('x');
% ylabel('y');

