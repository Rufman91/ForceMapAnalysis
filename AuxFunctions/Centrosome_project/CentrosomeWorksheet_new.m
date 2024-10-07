%% Tip data cleanup and deconvolution DO THIS WITH THE OTHER SCRIPT

Idx = 2;

InChannel = E.CantileverTips{Idx}.get_channel('Height (measured)');

OutChannel = AFMBaseClass.set_freehand_area_in_channel_to_value(InChannel, min(InChannel.Image,[],'all'));

OutChannel.Name = 'Processed';

E.CantileverTips{Idx}.add_channel(OutChannel,true)

E.CantileverTips{1}.deconvolute_cantilever_tip(1e-9,'Processed',false);
E.CantileverTips{2}.deconvolute_cantilever_tip(1e-9,'Processed',false);

%% Set nominal bead radius

for i=1:E.NumForceMaps
    E.FM{i}.TipRadius = 1e-6;
end

%% Fit force curves with simple models to get Contact Height

% Contact Point option: Fast
% SensitivityCorrectionMethod: KeepOriginal
% AllowXShift: true

% E.force_map_analysis_general

for i=1:E.NumForceMaps
E.FM{i}.CP = E.FM{i}.CP_HertzFitted;
E.FM{i}.create_height_map_by_current_cp_and_level_by_other_channel('Height (measured)',false)
end

%% Centrosomes smoothing steps

for i =1:E.NumForceMaps
    Chan = E.FM{i}.get_channel('Contact Height');
    Map = Chan.Image;
    Map = AFMImage.find_and_replace_outlier_lines(Map,3);
    Chan.Image = Map;
    SmoothChan = Chan;
    SmoothChan.Image = sgolayfilt(SmoothChan.Image,3,15);
    figure
    imshowpair(Chan.Image,SmoothChan.Image,'montage')
    SmoothChan.Name = 'Contact Height Smoothed';
    E.FM{i}.add_channel(SmoothChan,true);
end
%% Background selection

E.choose_segments_manually

%% Bias correction

for i=1:E.NumForceMaps
    Dat = E.FM{i}.get_segment_data_from_channel('Contact Height Smoothed', 'MatchString', 'Seg-01');
    Bias = mean(Dat(:,1),'all','omitnan'); 
    Chan = E.FM{i}.get_channel('Contact Height Smoothed');
    Chan.Image = Chan.Image - Bias;
    E.FM{i}.add_channel(Chan,true);
end

%% Now for thin film analysis
% Choose a thin film model in set_force_map_analysis_options
% Leave everything else as is

E.set_force_map_analysis_options

E.ForceMapAnalysisOptions.EModOption.Hertz.TopographyHeightChannel = 'Contact Height Smoothed';

% Analyze!
E.force_map_analysis_general


%% Create an artificial AFM tip for devonvolution

% E.create_artificial_cantilever_tip('HalfsphereTip','halfsphere','Radius',1e-6,'TipHeight',2e-6,'ImageResolution',128)

E.create_artificial_cantilever_tip('HalfsphereTip','halfsphere',...
    'Radius',1e-6,'TipHeight',2e-6,'ImageResolution',128,...
    'FuseMethod','linear',...
    'TipApexChannel',[],...
    'HeightDifference',33e-9);
AFMImage.plot_mesh_channel_to_scale(E.CantileverTips{end}.get_channel('Eroded Tip'))
E.create_artificial_cantilever_tip('HalfsphereTip','halfsphere',...
    'Radius',1e-6,'TipHeight',2e-6,'ImageResolution',128,...
    'FuseMethod','linear',...
    'TipApexChannel',E.CantileverTips{2}.get_channel('Eroded Tip'),...
    'HeightDifference',33e-9);
AFMImage.plot_mesh_channel_to_scale(E.CantileverTips{end}.get_channel('Eroded Tip'))

% %% Get Topographical data an metrics
% 
% for i=1:E.NumForceMaps
%     i
%     if i < 9
%         CTIndex = 1;
%     else
%         CTIndex = 2;
%     end
%     R = 1e-6;
%     Radii = E.FM{i}.convert_data_list_to_map(min(sqrt(2*R*E.FM{i}.IndentationDepthHertz - E.FM{i}.IndentationDepthHertz.^2),R));
%     
%     E.FM{i}.localSurfaceFit_ClassWrapper('Deconvoluted',Radii,false)
%     
% end
% 
% %% Deconvoluting Images
% 
% for i=1:E.NumForceMaps
%     i
%     if i <= 8
%         TipIdx = 6;
%     else
%         TipIdx = 6;
%     end
% %     TipIdx = 4;
%     E.FM{i}.deconvolute_image(E.CantileverTips{TipIdx},'Contact Height Smoothed',1024,true);  
% end

%% Indentation modulus readout after having segmented Seg 01 as Background 
% and Seg 02 as Centrosome. 'IncludeIndexVector',[2] makes sure, just
% centrosome data is read out
% Indentation modulus readout and calculations

% Initialize cell arrays to store results
Height = cell(1, E.NumForceMaps);
Slope = cell(1, E.NumForceMaps);
EffIndRad = cell(1, E.NumForceMaps);
EModSimple = cell(1, E.NumForceMaps);
EModK1000 = cell(1, E.NumForceMaps);
EModK250 = cell(1, E.NumForceMaps);
EModThinFilm = cell(1, E.NumForceMaps);
Volumes = zeros(1, E.NumForceMaps);
EquivalentRadii = zeros(1, E.NumForceMaps);
FilteredData = struct('EModSimple', {}, 'EModK1000', {}, 'EModK250', {}, 'EModThinFilm', {}, 'Height', {}, 'Volume', {});

% Loop through each Force Map
for i = 1:E.NumForceMaps
    Height{i} = E.FM{i}.get_segment_data_from_channel('Contact Height Smoothed', 'MatchString', 'Seg-02');
%     Slope{i} = E.FM{i}.get_segment_data_from_channel('Local Slope Kernel-R. = 2.50e-07 (01)', 'MatchString', 'Seg-02');
    EffIndRad{i} = E.FM{i}.get_segment_data_from_channel('Effective Radius Hertz (02)', 'MatchString', 'Seg-02');
    EModSimple{i} = E.FM{i}.get_segment_data_from_channel('Indentation Modulus Hertz', 'MatchString', 'Seg-02');
    EModThinFilm{i} = E.FM{i}.get_segment_data_from_channel('Indentation Modulus Hertz (02)', 'MatchString', 'Seg-02');
%     EModK1000{i} = real(E.FM{i}.get_segment_data_from_channel('Indentation Modulus Hertz (07)', 'MatchString', 'Seg-02'));
%     EModK250{i} = real(E.FM{i}.get_segment_data_from_channel('Indentation Modulus Hertz (09)', 'MatchString', 'Seg-02'));

    % Calculate volume (treat negative heights as zero)
    positiveHeight = max(Height{i}, 0);
    Volume = sum(positiveHeight) * (E.FM{i}.ScanSizeX/E.FM{i}.NumPixelsX * E.FM{i}.ScanSizeY/E.FM{i}.NumPixelsY);
    Volumes(i) = Volume;

    % Equivalent radius of a sphere with the same volume
    EquivalentRadii(i) = (3 * Volume / (4 * pi))^(1/3);

    % Filter data
    maxHeight = max(Height{i});
    validIndices = Slope{i} <= 0.25 & Height{i} > 0.8 * maxHeight & ...
                   ~isnan(Height{i}) & ~isnan(Slope{i}) & ~isnan(EffIndRad{i}) & ...
                   ~isnan(EModSimple{i}) & ~isnan(EModK1000{i}) & ~isnan(EModK250{i}) & ~isnan(EModThinFilm{i});
    
    FilteredData(i).EModSimple = EModSimple{i}(validIndices);
%     FilteredData(i).EModK1000 = EModK1000{i}(validIndices);
%     FilteredData(i).EModK250 = EModK250{i}(validIndices);
    FilteredData(i).EModThinFilm = EModThinFilm{i}(validIndices);
    FilteredData(i).Height = Height{i}(validIndices);
    FilteredData(i).Volume = EquivalentRadii(i);
end


% Calculate means and standard deviations
meanEModSimple = arrayfun(@(i) mean(FilteredData(i).EModSimple), 1:E.NumForceMaps);
stdEModSimple = arrayfun(@(i) std(FilteredData(i).EModSimple), 1:E.NumForceMaps);
meanEModThinFilm = arrayfun(@(i) mean(FilteredData(i).EModThinFilm), 1:E.NumForceMaps);
stdEModThinFilm = arrayfun(@(i) std(FilteredData(i).EModThinFilm), 1:E.NumForceMaps);
meanEModK1000 = arrayfun(@(i) mean(FilteredData(i).EModK1000), 1:E.NumForceMaps);
stdEModK1000 = arrayfun(@(i) std(FilteredData(i).EModK1000), 1:E.NumForceMaps);
meanEModK250 = arrayfun(@(i) mean(FilteredData(i).EModK250), 1:E.NumForceMaps);
stdEModK250 = arrayfun(@(i) std(FilteredData(i).EModK250), 1:E.NumForceMaps);

RobustFit = false;

% 2. Mean EMod vs Centrosome Volume
figure;
subplot(4, 1, 1);
errorbar(EquivalentRadii, meanEModSimple, stdEModSimple, 'o');
xlabel('Centrosome Equivalent Radius');
ylabel('Mean EMod Simple');
title('Mean EMod Simple vs Centrosome Volume');
addCorrelationInfo(EquivalentRadii, meanEModSimple, RobustFit);

subplot(4, 1, 2);
errorbar(EquivalentRadii, meanEModThinFilm, stdEModThinFilm, 'o');
xlabel('Centrosome Equivalent Radius');
ylabel('Mean EMod ThinFilm');
title('Mean EMod ThinFilm vs Centrosome Volume');
addCorrelationInfo(EquivalentRadii, meanEModThinFilm, RobustFit);

subplot(4, 1, 3);
errorbar(EquivalentRadii, meanEModK1000, stdEModK1000, 'o');
xlabel('Centrosome Equivalent Radius');
ylabel('Mean EMod K1000');
title('Mean EMod K1000 vs Centrosome Volume');
addCorrelationInfo(EquivalentRadii, meanEModK1000, RobustFit);

subplot(4, 1, 4);
errorbar(EquivalentRadii, meanEModK250, stdEModK250, 'o');
xlabel('Centrosome Equivalent Radius');
ylabel('Mean EMod K250');
title('Mean EMod K250 vs Centrosome Volume');
addCorrelationInfo(EquivalentRadii, meanEModK250, RobustFit);

% 3. Mean EMod vs Centrosome Height
maxHeight = arrayfun(@(i) max(FilteredData(i).Height), 1:E.NumForceMaps);

figure;
subplot(4, 1, 1);
errorbar(maxHeight, meanEModSimple, stdEModSimple, 'o');
xlabel('Centrosome Maximum Height');
ylabel('Mean EMod Simple');
title('Mean EMod Simple vs Centrosome Height');
addCorrelationInfo(maxHeight, meanEModSimple, RobustFit);

subplot(4, 1, 2);
errorbar(maxHeight, meanEModThinFilm, stdEModThinFilm, 'o');
xlabel('Centrosome Maximum Height');
ylabel('Mean EMod ThinFilm');
title('Mean EMod ThinFilm vs Centrosome Height');
addCorrelationInfo(maxHeight, meanEModThinFilm, RobustFit);

subplot(4, 1, 3);
errorbar(maxHeight, meanEModK1000, stdEModK1000, 'o');
xlabel('Centrosome Maximum Height');
ylabel('Mean EMod K1000');
title('Mean EMod K1000 vs Centrosome Height');
addCorrelationInfo(maxHeight, meanEModK1000, RobustFit);

subplot(4, 1, 4);
errorbar(maxHeight, meanEModK250, stdEModK250, 'o');
xlabel('Centrosome Maximum Height');
ylabel('Mean EMod K250');
title('Mean EMod K250 vs Centrosome Height');
addCorrelationInfo(maxHeight, meanEModK250, RobustFit);


function addCorrelationInfo(x, y, useRobustFit)
    % Calculate correlation coefficient and p-value
    [r, p] = corrcoef(x, y);
    r_value = r(1, 2);
    p_value = p(1, 2);
    
    % Calculate the number of samples needed for 80% power at alpha = 0.05
    Z_alpha_over_2 = 1.96; % for alpha = 0.05
    Z_beta = 0.84; % for 80% power
    fisher_z = 0.5 * log((1 + r_value) / (1 - r_value));
    needed_samples = (Z_alpha_over_2 + Z_beta)^2 / fisher_z^2 + 3;
    
    % Display correlation line
    hold on;
    if useRobustFit
        coeffs = robustfit(x, y);
        fittedX = linspace(min(x), max(x), 200);
        fittedY = coeffs(1) + coeffs(2) * fittedX;
    else
        coeffs = polyfit(x, y, 1);
        fittedX = linspace(min(x), max(x), 200);
        fittedY = polyval(coeffs, fittedX);
    end
    plot(fittedX, fittedY, 'r-', 'LineWidth', 1);
    
    % Display r, p values, and needed samples
    text(min(x) + 0.05 * (max(x) - min(x)), max(y) - 0.05 * (max(y) - min(y)), ...
         sprintf('r = %.2f\np = %.2g\nNeeded samples = %.0f', r_value, p_value, ceil(needed_samples)), ...
         'EdgeColor', 'black', 'BackgroundColor', 'white', 'Margin', 5);
    hold off;
end
