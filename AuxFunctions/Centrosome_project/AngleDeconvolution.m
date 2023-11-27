% centrosome AFM data analysis: angle deconvolution to identify areas from
% which indentation modulus values are extractedE.sh% post-processing step after having evaluated the quality of the fits and
% sensitivity correction method using E.show_image
% Julia Garcia Baucells 2022

clear
E = Experiment.load;
cd(E.ExperimentFolder)

dbstop in AngleDeconvolution.m at 198
format long g;
format compact;
workspace;  % make sure the workspace panel is showing

close all
clc
show_fig = 'on';
msg1 = "Do you need to exclude force maps?";
opts1 = ["Yes" "No"];
choice1 = menu(msg1,opts1);

msg2 = "Do you want to apply a subsequent ForceMapAnalysisOptions?";
opts2 = ["Yes" "No"];
choice2 = menu(msg2,opts2);

if choice2 == 1
    msg3 = "Which ForceMapAnalysisOptions do you want to apply?";
    opts3 = ["01" "02" "03"];
    choice3 = menu(msg3,opts3);

    s1 = ' Fit Range';
    s2 = ' ('+opts2(choice3)+')';
else
    s1 = '';
    s2 = '';
    choice3 = []; 
end

% read force map number to be excluded
if choice1 == 1
    if isfile('Exclude.txt')
        fileID = fopen('Exclude.txt', 'r');
        datacell = textscan(fileID, '%f', 'Delimiter',' ', 'CollectOutput', 1);
        fclose(fileID);
        datavalues = unique(datacell{1}); % Remove duplicates
    else
        datavalues = [];
    end
    if choice3 == 3
        if isfile('Exclude03.txt')
            fileID = fopen('Exclude03.txt', 'r');
            datacell03 = textscan(fileID, '%f', 'Delimiter',' ', 'CollectOutput', 1);
            fclose(fileID);
            datavalues03 = unique(datacell03{1}); % Remove duplicates

            % combine values without repeating
            datavalues = unique([datavalues; datavalues03]);
        end
    end
else
    datavalues = [];
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
                [PlaneAngle] = calculate_angle_topography(arr);
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
        ChannelEModHertz = E.FM{m}.get_channel(strcat('Indentation Modulus Hertz',s2));
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
        ChannelContact = E.FM{m}.get_channel(strcat('Contact Height',s2));
        CsFlatHeight = ChannelContact.Image.*AngleCsBW.*QFits;
        CsFlatHeight(CsFlatHeight==0) = NaN;
        
        ChannelIndenDepth = E.FM{m}.get_channel(strcat('Indentation Depth Hertz',s1,s2));
        CsFlatInden = ChannelIndenDepth.Image.*AngleCsBW.*QFits;
        CsFlatInden(CsFlatInden==0) = NaN;

        % volume approximation
        CsContactHeight = ChannelContact.Image.*erodedBW2; 
        CsVolume_voxel = CsContactHeight.*pxSize^2; 
        CsVolume = sum(CsVolume_voxel(:)); 

        % save 
        cd(E.ForceMapFolders{m,1})
        filename = strcat('Processed',s2);
        save(filename,'T2','CsArea','CsRadius','CsEModHertz','CsFlatArea','CsFlatHeight','CsFlatInden','CsVolume')

        cd ..
        close all
    end
end