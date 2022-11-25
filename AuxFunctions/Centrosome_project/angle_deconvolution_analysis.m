% centrosome AFM data analysis: angle deconvolution to identify areas from
% which indentation modulus values are extracted
% post-processing step after having evaluated the quality of the fits and
% sensitivity correction method using E.show_image
% Julia Garcia Baucells 2022

% clear
% E = Experiment.load;
cd(E.ExperimentFolder)

format long g;
format compact;
workspace;  % make sure the workspace panel is showing

close all
clc
show_fig = 'on';

% read force map number to be excluded
if isfile('Exclude.txt')
    fileID = fopen('Exclude.txt', 'r');
    datacell = textscan(fileID, '%f', 'Delimiter',' ', 'CollectOutput', 1);
    fclose(fileID);
    datavalues = datacell{1};    %as a numeric array
end

for m = 1:E.NumForceMaps
    if ismember(m, datavalues) 
        % skip evaluation round
    else
        % height data segmentation
        ChannelHeight = E.FM{m}.get_channel('Processed');
        figure('name','Processed','visible',show_fig); hold on
        imagesc(ChannelHeight.Image); axis image; colorbar
        T = multithresh(ChannelHeight.Image);
        BW = imbinarize(ChannelHeight.Image,T);
        BW2 = imfill(BW,'holes');
        SE = strel("disk",2);
        erodedBW2 = imerode(BW2,SE);
        dilatedBW2 = imdilate(BW2,SE);
        figure('name','Processed segmented','visible',show_fig); hold on
        imagesc(erodedBW2); axis image

        % get indentation depth average
        ChannelInden = E.FM{m}.get_channel('Indentation Depth Hertz');
        CsInden = ChannelInden.Image.*erodedBW2;
        figure('name', 'Indentation Depth','visible',show_fig); hold on
        imagesc(CsInden); axis image; colorbar
        CsInden(CsInden==0)=NaN;
        CsIndenAvrg = mean(CsInden(:), 'omitnan');

        % tip radius at this identation depth
        TipAreaFX = E.CantileverTips{1}.ProjectedTipArea;
        TipAreaCsInden = TipAreaFX(round(CsIndenAvrg*1e+9));
        TipRadiusCsInden = sqrt(TipAreaCsInden/pi);

        % px size
        pxSize = (ChannelHeight.ScanSizeX/ChannelHeight.NumPixelsX);

        CsProps = regionprops(erodedBW2, 'Area','EquivDiameter');
        [~, ind] = sort([CsProps.Area], 'descend');
        CsProps = CsProps(ind);
        CsArea = CsProps(1).Area.*pxSize^2;
        CsRadius = (CsProps(1).EquivDiameter*pxSize)/2;

        % calculate size of struct element
        SeRadius = round(TipRadiusCsInden/pxSize);
        kernel = strel('disk', SeRadius);
        dimensions = size(ChannelHeight.Image);
        dimensions2 = size(kernel.Neighborhood);

        % define kernel center indices
        kernelCenter_x = round(dimensions2(1)/2);
        kernelCenter_y = round(dimensions2(2)/2);

        AngleImage = zeros(dimensions(1),dimensions(2));
        for i = 1:dimensions(1)
            for j = 1:dimensions(2)
                for k = 1:dimensions2(1)
                    for l = 1:dimensions2(2)
                        ii = i+(k-kernelCenter_x);
                        jj = j+(l-kernelCenter_y);
                        if (ii >= 1 && ii <= dimensions(1) && jj >= 1 && jj <= dimensions(2))
                            Cloud_DPs(ii,jj) = ChannelHeight.Image(ii,jj)*kernel.Neighborhood(k,l); % data points within kernel
                            %                         AngleImage(i,j) = AngleImage(i,j) + ChannelHeight.Image(ii,jj)*kernel.Neighborhood(k,l);
                        end
                    end
                end
                %         imagesc(zeros(dimensions(1),dimensions(2))); axis image
                %         hold on
                %         imagesc(Cloud_DPs)
                Cloud_DPs(Cloud_DPs == 0) = NaN;
                [row, col] = ind2sub(size(Cloud_DPs), 1:numel(Cloud_DPs));
                arr = [row(:).*pxSize, col(:).*pxSize, Cloud_DPs(:)]; % convert axis units to m
                arr(any(isnan(arr), 2), :) = []; % eliminate rows with nans
                [PlaneAngle] = calculate_topography_angle(arr);
                AngleImage(i,j) = PlaneAngle;
                clear Cloud_DPs
            end
        end
        figure('name', 'Topography angle')
        imagesc(flipud(AngleImage)); axis image; colorbar

        % angle-based segmentation
        T2 = multithresh(AngleImage);
        AngleCsBW = imbinarize(AngleImage,T2).*erodedBW2; % exclude flat areas outside centrosome
        AngleBkgBW = imbinarize(AngleImage,T2).*imcomplement(dilatedBW2); % glass background area
        figure('name', 'Topography angle segmented','visible',show_fig)
        subplot(1,2,1)
        hold on
        imagesc(AngleCsBW); axis image; title('centrosome')
        subplot(1,2,2)
        hold on
        imagesc(AngleBkgBW); axis image; title('glass')

        % get Indentation Modulus Hertz data
        ChannelEModHertz = E.FM{m}.get_channel('Indentation Modulus Hertz');
        % centrosome
        CsEModHertz = ChannelEModHertz.Image.*AngleCsBW;
        CsEModHertz(CsEModHertz==0)=NaN;
        figure('name', 'Indentation Modulus Hertz')
        subplot(1,2,1)
        hold on
        imagesc(CsEModHertz*1e-3); axis image; colorbar; title('centrosome') % kPa
        % glass
        BkgEModHertz = ChannelEModHertz.Image.*AngleBkgBW;
        BkgEModHertz(BkgEModHertz==0)=NaN;
        subplot(1,2,2)
        hold on
        imagesc(BkgEModHertz*1e-3); axis image; colorbar; title('glass') % kPa

        % calculate flat area
        CsFlatProps = regionprops(AngleCsBW, 'Area');
        CsFlatArea = CsFlatProps.Area*pxSize^2;
        BkgFlat = regionprops(AngleBkgBW, 'Area');
        BkgFlatArea = BkgFlat.Area*pxSize^2;

        % centrosome height + indentation on the flat area
%         ChannelContact = E.FM{m}.get_channel('Contact Height');
        CsFlatHeight = ChannelHeight.Image.*AngleCsBW;
        CsFlatHeight(CsFlatHeight==0)=NaN;
        CsFlatInden = ChannelInden.Image.*AngleCsBW;
        CsFlatInden(CsFlatInden==0)=NaN;

        cd(E.ForceMapFolders{m,1})
        filename = 'Processed';
        save(filename,'T2','CsEModHertz','BkgEModHertz','CsArea','CsRadius','CsFlatArea','BkgFlatArea','CsFlatHeight',...
        'CsFlatInden')

        cd ..
        close all
    end
end