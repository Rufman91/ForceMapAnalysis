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
msg = "Do you need to exclude force maps?";
opts = ["Yes" "No"];
choice = menu(msg,opts);

% read force map number to be excluded
if choice == 1
    if isfile('Exclude.txt')
        fileID = fopen('Exclude.txt', 'r');
        datacell = textscan(fileID, '%f', 'Delimiter',' ', 'CollectOutput', 1);
        fclose(fileID);
        datavalues = datacell{1};    %as a numeric array
    end
else 
    datavalues = 0; 
end

for m = 1:E.NumForceMaps
    if ismember(m, datavalues) 
        % skip evaluation round
    else
        % processed height data segmentation (centrosome 2D-stamp on glass) 
        ChannelProcHeight = E.FM{m}.get_channel('Processed');
        figure('name','Processed','visible',show_fig); hold on
        imagesc(ChannelProcHeight.Image.*1e9); axis image; 
        c = colorbar; c.Location = 'northoutside'; c.Label.String = 'Processed [nm]'; 
        set(gca,'FontSize', 16, 'Linewidth', 1.5); axis off
        T = multithresh(ChannelProcHeight.Image);
        BW = imbinarize(ChannelProcHeight.Image,T);
        BW2 = imfill(BW,'holes');
        SE = strel("disk",2);
        erodedBW2 = imerode(BW2,SE);
        dilatedBW2 = imdilate(BW2,SE);
        figure('name','Processed mask','visible',show_fig); hold on
        imagesc(erodedBW2); axis image; axis off 

        % filter out bad R2 Hertz fits 
        PredictiveR2Thr = 0.96; 
        R2Thr = 0.99; 
        ChannelPredictiveR2 = E.FM{m}.get_channel('Hertz Fit Predictive RSquare');
        ChannelR2 = E.FM{m}.get_channel('Hertz Fit RSquare');
        QPredictiveR2 = ChannelPredictiveR2.Image; 
        QPredictiveR2(QPredictiveR2<PredictiveR2Thr) = 0; 
        QR2 = ChannelR2.Image; 
        QR2(QR2<R2Thr) = 0; QFits = QR2.*QPredictiveR2; QFits(QFits>0)=1; 

        % get indentation depth average in centrosome region
        ChannelIndenDepth = E.FM{m}.get_channel('Indentation Depth Hertz');
        CsInden = ChannelIndenDepth.Image.*erodedBW2.*QFits;
        figure('name', 'Indentation depth','visible',show_fig); hold on
        imagesc(CsInden.*1e9); axis image; c = colorbar; 
        c.Location = 'northoutside'; c.Label.String = 'Indentation Depth [nm]'; 
        set(gca,'FontSize', 16, 'Linewidth', 1.5); axis off
        CsInden(CsInden==0)=NaN; CsIndenAvrg = mean(CsInden(:), 'omitnan'); 

        % tip radius at this identation depth
        TipAreaFX = E.CantileverTips{1}.ProjectedTipArea;
        TipAreaCsInden = TipAreaFX(round(CsIndenAvrg*1e+9));
        TipRadiusCsInden = sqrt(TipAreaCsInden/pi);

        % px size
        pxSize = (ChannelProcHeight.ScanSizeX/ChannelProcHeight.NumPixelsX);

        CsProps = regionprops(erodedBW2, 'Area','EquivDiameter');
        [~, ind] = sort([CsProps.Area], 'descend');
        CsProps = CsProps(ind);
        CsArea = CsProps(1).Area.*pxSize^2;
        CsRadius = (CsProps(1).EquivDiameter*pxSize)/2;

        % calculate size of struct element
        SeRadius = round(TipRadiusCsInden/pxSize);
        kernel = strel('disk', SeRadius);
        dimensions = size(ChannelProcHeight.Image);
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
                            Cloud_DPs(ii,jj) = ChannelProcHeight.Image(ii,jj)*kernel.Neighborhood(k,l); % data points within kernel
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
        figure('name', 'Topography angle', 'visible', show_fig);
        imagesc(flipud(AngleImage)); axis image; c = colorbar; 
        c.Location = 'northoutside'; c.Label.String = 'Angle [rad]'; % 90 deg flat 
        set(gca,'FontSize', 16, 'Linewidth', 1.5); axis off
        
        % angle-based segmentation
        T2 = multithresh(AngleImage);
        AngleCsBW = imbinarize(AngleImage,T2).*erodedBW2; % exclude flat areas outside centrosome
        AngleBkgBW = imbinarize(AngleImage,T2).*imcomplement(dilatedBW2); % glass background area
        figure('name', 'Topography angle segmented','visible',show_fig)
        subplot(1,2,1); hold on
        imagesc(AngleCsBW); axis image; title('centrosome'); axis off 
        set(gca,'FontSize', 16, 'Linewidth', 1.5); 
        subplot(1,2,2); hold on
        imagesc(AngleBkgBW); axis image; title('glass'); axis off 
        set(gca,'FontSize', 16, 'Linewidth', 1.5); 

        % get indentation modulus Hertz data 
        ChannelEModHertz = E.FM{m}.get_channel('Indentation Modulus Hertz');
        % from centrosome's flat region
        CsEModHertz = ChannelEModHertz.Image.*AngleCsBW.*QFits;
        CsEModHertz(CsEModHertz==0)=NaN;
        figure('name', 'Indentation modulus Hertz')
        hold on
        imagesc(CsEModHertz*1e-3); axis image; c = colorbar; % kPa
        c.Location = 'northoutside'; c.Label.String = 'Indentation modulus Hertz [kPa]'; 
        axis off; set(gca,'FontSize', 16, 'Linewidth', 1.5); 

        % calculate centrosome's flat area
        CsFlatProps = regionprops(AngleCsBW, 'Area');
        CsFlatArea = CsFlatProps.Area*pxSize^2;

        % centrosome height (contact height) & indentation on the flat area
        ChannelContact = E.FM{m}.get_channel('Contact Height');
        CsFlatHeight = ChannelContact.Image.*AngleCsBW.*QFits;
        CsFlatHeight(CsFlatHeight==0) = NaN;
        
        CsFlatInden = ChannelIndenDepth.Image.*AngleCsBW.*QFits;
        CsFlatInden(CsFlatInden==0) = NaN;

        % save 
        cd(E.ForceMapFolders{m,1})
        filename = 'Processed';
        save(filename,'T2','CsArea','CsRadius','CsEModHertz','CsFlatArea','CsFlatHeight','CsFlatInden')

        cd ..
        close all
    end
end