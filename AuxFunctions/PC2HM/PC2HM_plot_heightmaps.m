function PC2HM_plot_heightmaps(string, page)
if nargin < 2
    page = 1;
end
% Get all files in the current folder that contain the input string
files = dir(['*', string, '*']);

% Determine number of files and choose tiledlayout dimensions
num_files = length(files);
pages = ceil(num_files/4);

% Create a figure
figure('units','normalized','outerposition',[0 0 1 1])

i = page;
% Create a tiledlayout with the determined dimensions
t = tiledlayout(2, 2);
for j = 1:4
    idx = (i-1)*4 + j;
    if idx > num_files
        break;
    end
    nexttile;
    clear HeightMap Lambda Sigma Noise MaxPeakMap MinPeakMap MaxPeakValueMap;
    load(files(idx).name,'HeightMap', 'MaxPeakMap', 'MinPeakMap', 'MaxPeakValueMap','Lambda', 'Sigma', 'Noise');
    if ~exist('HeightMap')
        HeightMap = MaxPeakMap;
    end
    I = imshow(HeightMap,[]);
    colormap(AFMImage.define_afm_color_map);
    lambda = Lambda;
    sigma = Sigma;
    noise = Noise;
    title({files(idx).name, ['Lambda: ', num2str(lambda)], ['Sigma: ', num2str(sigma)], ['Noise: ', num2str(noise)]});
end
if i < pages
    uicontrol('Style', 'pushbutton', 'String', 'Next',...
        'Units','normalized','Position', [0.9 0.1 0.1 0.05],...
        'Callback', {@next_page, i+1});
end
if i > 1
    uicontrol('Style', 'pushbutton', 'String', 'Back',...
        'Units','normalized','Position', [0.8 0.1 0.1 0.05],...
        'Callback', {@prev_page, i-1});
end

    function next_page(src, event, page)
        close(gcf);
        PC2HM_plot_heightmaps(string, page);
    end
    function prev_page(src, event, page)
        close(gcf);
        PC2HM_plot_heightmaps(string, page);
    end

end

