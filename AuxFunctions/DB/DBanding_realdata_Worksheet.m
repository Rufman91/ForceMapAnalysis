clear Skew noisy_signal
clear peak_skewness_minus peak_skewness_plus period amplitude y_filtered

lower_cutoff = 10e-9;
upper_cutoff = 80e-9;
padding_factor = 2;

NI = 3;
NSegs = numel(E.I{NI}.Segment);

for i=1:NSegs/2
    X = E.I{NI}.Segment(NSegs/2 + i).RelativePosition;
    Y = E.I{NI}.Segment(NSegs/2 + i).Height;
    
    noisy_signal_final = detrend(Y,5);
    XQ = X;
%     XQ = linspace(min(X),max(X),1000);
%     noisy_signal_final = interp1(X,noisy_signal_final,XQ,'makima');
    [period(i), amplitude(i), y_filtered] = find_spatial_frequency(XQ, noisy_signal_final, lower_cutoff, upper_cutoff, padding_factor,true);
    peak_skewness_plus = periodic_signal_skewness(XQ, y_filtered, period(i));
    peak_skewness_minus = -periodic_signal_skewness(XQ, -y_filtered, period(i));
    
    Skew(i) = median([peak_skewness_plus ; peak_skewness_minus]);
end
figure
subplot(3,1,1)
plot(XQ,noisy_signal_final)
hold on
legend
subplot(3,1,2)
plot(period)
title(['Median period = ' num2str(median(period))])
xticks([1:NSegs/2+1])
xticklabels({E.I{NI}.Segment(NSegs/2+1:end).Name})
subplot(3,1,3)
plot(Skew','-O')
legend
hold on
plot((range(Skew)+.1)/2*sign(Skew'),'-O')
xticks([1:NSegs/2+1])
xticklabels({E.I{NI}.Segment(NSegs/2+1:end).Name})