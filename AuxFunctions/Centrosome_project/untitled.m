for i=1:E.NumForceMaps
    Dat = E.FM{i}.get_segment_data_from_channel('Contact Height Smoothed');
    Bias = mean(Dat,'all','omitnan'); 
    Chan = E.FM{i}.get_channel('Contact Height Smoothed');
    Chan.Image = Chan.Image -Bias;
    E.FM{i}.add_channel(Chan,true);
end
