% close all 
% clc 

ChanHeight = E.FM{1}.get_channel('Processed'); 
figure(); hold on 
imagesc(ChanHeight.Image); axis image
SE = strel('disk', 5);
FiltChanHeight = stdfilt(ChanHeight.Image, SE.Neighborhood); % standard deviation filtering, 
% the value of each output pixel is the standard deviation of the 5-by-5 neighborhood (disk) 
% around the corresponding input pixel - high signal areas have high variability
figure(); hold on 
imagesc(FiltChanHeight); axis image; colorbar
% Segmentation based on region homogeneity 
AreasChanHeight = FiltChanHeight; 
Thr = 1.2e-8; 
AreasChanHeight(AreasChanHeight > Thr) = 0;
AreasChanHeight(AreasChanHeight>0) = 1; % binarize 
figure(); hold on 
imagesc(AreasChanHeight); axis image;

% Segmentation of the processed height data - get only areas inside
% centrosome
T = multithresh(ChanHeight.Image);
Cs_Seg = imbinarize(ChanHeight.Image,T); % works well 
figure(); hold on 
imagesc(Cs_Seg); axis image 
AreasCs_Seg = AreasChanHeight.*Cs_Seg; 
figure(); hold on 
imagesc(AreasCs_Seg); axis image 
