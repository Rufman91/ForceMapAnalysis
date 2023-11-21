format long g;
format compact;
workspace;  % make sure the workspace panel is showing
% dbstop in collect_plot_results.m at 130
% dbstop in collect_plot_results.m at 175

close all
clear
clc

path1 = '/Users/julia/MyFolder/CentrosomesExchange/Mjtraining/';
addpath('internal functions');

DirOutput = dir(path1);
FileNames = {DirOutput.name}';
NumberFolders = numel(FileNames);

opts = ["", "01", "02"]; 
pat = 'colloidal';

% filter out unwanted folders
counter = 1;
sample_idx = nan(NumberFolders,1);
for m = 1:NumberFolders
    TF = contains(FileNames(m,1),pat);
    if TF == 1
        sample_idx(counter) = m;
    else
        sample_idx(counter) = 0;
    end
    counter = counter + 1;
end
sample_idx = nonzeros(sample_idx);
clear TF; clear counter;

counter(1:length(opts), 1:length(sample_idx)) = 0; 
for s = 1:length(opts)
for n = 1:length(sample_idx)
    ii = sample_idx(n);
    path2 = char(strcat(path1,'/',cellstr(FileNames(ii,1)),'/','ZoomSweepFast'));
    cd(path2); % go inside date folder

    load('ZoomSweepFast.mat')
    if isfile('Exclude.txt')
        fileID = fopen('Exclude.txt', 'r');
        datacell = textscan(fileID, '%f', 'Delimiter',' ', 'CollectOutput', 1);
        fclose(fileID);
        datavalues = datacell{1};    %as a numeric array
    end
   
    % collect data
    for i = 1:obj.NumForceMaps
        if ismember(i, datavalues)
            % skip evaluation round
        else
            cd(obj.ForceMapFolders{i,1})
            if isempty(isspace(opts(s)))
                s2 = ''; 
            else 
                s2 = ' ('+opts(s)+')'; 
            end 
            load(strcat('Processed',s2,'.mat'))
            counter(s,n) = counter(s,n) + 1;
            CsEModHertz_data{s,n,counter(s,n)} = CsEModHertz(:);
            CsEModHertz_mean{s,n,counter(s,n)} = mean(CsEModHertz(:),'omitnan');
            CsHeight_mean{s,n,counter(s,n)} = mean(CsFlatHeight(:),'omitnan');
            CsInden_mean{s,n,counter(s,n)} = mean(CsFlatInden(:),'omitnan');
            CsRadius_data{s,n,counter(s,n)} = CsRadius;
            CsAspectRatio{s,n,counter(s,n)} = mean(CsFlatHeight(:),'omitnan')/(CsRadius*2); 
            AngleThr{s,n,counter(s,n)} = T2;
            CsFlatArea_data{s,n,counter(s,n)} = CsFlatArea;
            CsVolume_data{s,n,counter(s,n)} = CsVolume; 
        end
    end
end
end

total_cs = sum(counter(1,:));
CsEModHertz_mean_l = nan(length(opts),total_cs);
Compression_l = nan(1,total_cs);
CsHeight_mean_l = nan(1,total_cs);
for s = 1:length(opts)
    for i = 1:length(sample_idx)
        for j = 1:counter(s,i)
            clr_counter = j + sum(counter(s,1:i-1));
            CsEModHertz_mean_l(s, clr_counter) = CsEModHertz_mean{s,i,j}.*1e-3;
            Compression_l(clr_counter) = (CsInden_mean{1,i,j}./CsHeight_mean{1,i,j})*100;
            CsHeight_mean_l(clr_counter) = CsHeight_mean{1,i,j}; 
        end
    end
end

figure(); hold on; 
CsEModHertz_std = std(CsEModHertz_mean_l,1)./mean(CsEModHertz_mean_l,1,'omitnan') ;
scatter(CsEModHertz_std, Compression_l, 70, 'b', "filled");
xlabel('Indentation modulus DV [%]');
ylabel('Compression [%]'); ylim([0 60]); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);

figure(); hold on; 
scatter(CsEModHertz_std, CsHeight_mean_l*1e9, 70, 'b', "filled");
xlabel('Indentation modulus DV [%]');
ylabel('Centrosome height [nm]'); hold on
box on; set(gca,'FontSize', 16, 'Linewidth', 1.5);
ylim([0 1000])

