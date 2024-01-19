%% Tip data cleanup and deconvolution DO THIS WITH THE OTHER SCRIPT
% InChannel = E.CantileverTips{1}.get_channel('Processed');
% for i=1:4
%     OutChannel = set_freehand_area_in_channel_to_value(InChannel, min(InChannel.Image,[],'all'));
% end
% 
% E.CantileverTips{1}.deconvolute_cantilever_tip;

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
i
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
    Dat = E.FM{i}.get_segment_data_from_channel('Contact Height Smoothed');
    Bias = mean(Dat,'all','omitnan'); 
    Chan = E.FM{i}.get_channel('Contact Height Smoothed');
    Chan.Image = Chan.Image-Bias;
    E.FM{i}.add_channel(Chan,true);
end

%% Now for thin film ana.
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

%% Get Topographical data an metrics
for i=1:E.NumForceMaps
    i
    if i < 9
        CTIndex = 1;
    else
        CTIndex = 2;
    end
    R = 1e-6;
    Radii = min(sqrt(2*R*E.FM{i}.IndentationDepthHertz - E.FM{i}.IndentationDepthHertz.^2),R);
    
    E.FM{i}.localSurfaceFit_ClassWrapper('Contact Height Smoothed',Radii,false)
    
end



%% Deconvoluting Images

for i=1:E.NumForceMaps
    i
    if i <= 8
        TipIdx = 6;
    else
        TipIdx = 6;
    end
%     TipIdx = 4;
    E.FM{i}.deconvolute_image(E.CantileverTips{TipIdx},'Contact Height Smoothed',1024,true);
    
end