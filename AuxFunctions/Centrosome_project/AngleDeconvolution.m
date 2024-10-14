% clear
% E = Experiment.load;
cd(E.ExperimentFolder)

dbstop in AngleDeconvolution.m at 111
format long g;
format compact;
workspace;  % Make sure the workspace panel is showing

close all
clc
show_fig = 'off';
msg1 = "Do you need to exclude force maps?";
opts1 = ["Yes" "No"];
choice1 = menu(msg1,opts1);

msg2 = "Do you want to apply a subsequent ForceMapAnalysisOptions?";
opts2 = ["Yes" "No"];
choice2 = menu(msg2,opts2);

if choice2 == 1
    msg3 = "Which ForceMapAnalysisOptions do you want to apply?";
    opts3 = ["01" "02" "03" "04" "05"];
    choice3 = menu(msg3,opts3);

    s1 = ' Fit Range';
    s2 = ' ('+opts3(choice3)+')';
else
    s1 = '';
    s2 = '';
    choice3 = []; 
end

% Read number of force maps to exclude
if choice1 == 1
    if isfile('Exclude.txt')
        fileID = fopen('Exclude.txt', 'r');
        datacell = textscan(fileID, '%f', 'Delimiter',' ', 'CollectOutput', 1);
        fclose(fileID);
        datavalues = unique(datacell{1}); % Remove duplicates
    else
        datavalues = [];
    end
%     if choice3 == 3
%         if isfile('Exclude03.txt')
%             fileID = fopen('Exclude03.txt', 'r');
%             datacell03 = textscan(fileID, '%f', 'Delimiter',' ', 'CollectOutput', 1);
%             fclose(fileID);
%             datavalues03 = unique(datacell03{1}); % Remove duplicates
% 
%             % combine values without repeating
%             datavalues = unique([datavalues; datavalues03]);
%         end
%     end
else
    datavalues = [];
end

for m = 1:E.NumForceMaps
    if ismember(m, datavalues) 
        % Skip evaluation round
    else
        % Contact height smoothed channel used for 2D centrosome segmentation  
        ChannelContactHeight = E.FM{m}.get_channel('Contact Height Smoothed');
        figure('name','Contact Height Smoothed','visible',show_fig); hold on
        imagesc(ChannelContactHeight.Image.*1e9); axis image; 
        c = colorbar; c.Location = 'northoutside'; c.Label.String = 'Height [nm]'; 
        set(gca,'FontSize', 16, 'Linewidth', 1.5); axis off; colormap("hot") 
        T = multithresh(ChannelContactHeight.Image);
        BW = imbinarize(ChannelContactHeight.Image,T);
        BW2 = imfill(BW,'holes');
        SE = strel("disk",2);
        erodedBW2 = imerode(BW2,SE);
        dilatedBW2 = imdilate(BW2,SE);
        figure('name','Centrosome mask','visible',show_fig); hold on
        imagesc(erodedBW2); axis image; axis off

        % Check if centrosome mask is open 
        logical_first_row = any(erodedBW2(1,:)); 
        logical_last_row = any(erodedBW2(end,:)); 
        logical_first_col = any(erodedBW2(:,1)); 
        logical_last_col = any(erodedBW2(:,end)); 
        folderPath = E.ForceMapFolders{m, 1};
        decisionFilePath = fullfile(folderPath, 'MaskUserDecision.mat');

        if any([logical_first_row, logical_last_row, logical_first_col, logical_last_col])
            if exist(decisionFilePath, 'file') == 2
                load(decisionFilePath, 'erodedBW2');
                figure('name','Centrosome mask','visible',show_fig); hold on
                imagesc(erodedBW2); axis image; axis off

            else
                userInput = input('Mask is open. Do you want to eliminate (1), keep (2), or trim manually (3)? ');
                switch userInput
                    case 1
                        % Eliminate the object at the image boundary
                        erodedBW2 = imclearborder(erodedBW2);
                        save(decisionFilePath, 'userInput', 'erodedBW2');
                    case 2
                        % Keep the object as it is
                         save(decisionFilePath, 'userInput', 'erodedBW2');
                    case 3
                        % Trim manually
                        roiwindow = CROIEditor_CE(erodedBW2);
                        waitfor(roiwindow, 'roi');
                        if ~isvalid(roiwindow)
                            disp('You closed ROI window without applying ROI')
                            return
                        end
                        [~,labels,~] = roiwindow.getROIData;
                        delete(roiwindow);
                        AppliedMask = (labels >= 1);
                        erodedBW2 = AppliedMask.*erodedBW2;
                        save(decisionFilePath, 'userInput', 'erodedBW2');
                end
            end
        else
            disp('Mask is closed. No need for user input.');
        end

        % Filter out bad R2 Hertz fits 
        PredictiveR2Thr = 0.96; 
        R2Thr = 0.99; 
        ChannelPredictiveR2 = E.FM{m}.get_channel('Hertz Fit Predictive RSquare');
        ChannelR2 = E.FM{m}.get_channel('Hertz Fit RSquare');
        QPredictiveR2 = ChannelPredictiveR2.Image; 
        QPredictiveR2(QPredictiveR2<PredictiveR2Thr) = 0; 
        QR2 = ChannelR2.Image; 
        QR2(QR2<R2Thr) = 0; QFits = QR2.*QPredictiveR2; QFits(QFits>0) = 1; 

        % Get indentation depth average in centrosome region
        ChannelIndenDepth = E.FM{m}.get_channel('Indentation Depth Hertz');
        CsInden = ChannelIndenDepth.Image.*erodedBW2.*QFits;
        figure('name', 'Indentation depth','visible',show_fig); hold on
        imagesc(CsInden.*1e9); axis image; c = colorbar; 
        c.Location = 'northoutside'; c.Label.String = 'Indentation Depth [nm]'; 
        set(gca,'FontSize', 16, 'Linewidth', 1.5); axis off
        CsInden(CsInden==0) = NaN; CsIndenAvrg = mean(CsInden(:), 'omitnan'); 

        % Tip radius at this identation depth
        if m>=1 && m<=45
            TipAreaFX = E.CantileverTips{1}.ProjectedTipArea;
        else
            TipAreaFX = E.CantileverTips{2}.ProjectedTipArea;
        end
        TipAreaCsInden = TipAreaFX(round(CsIndenAvrg*1e+9));
        TipRadiusCsInden = sqrt(TipAreaCsInden/pi);

        % Pixel size
        pxSize = (ChannelContactHeight.ScanSizeX/ChannelContactHeight.NumPixelsX);

        CsProps = regionprops(erodedBW2, 'Area','EquivDiameter');
        [~, ind] = sort([CsProps.Area], 'descend');
        CsProps = CsProps(ind);
        CsArea = CsProps(1).Area.*pxSize^2;
        CsRadiusXY = (CsProps(1).EquivDiameter*pxSize)/2;

        % Calculate size of struct element
        SeRadius = round(TipRadiusCsInden/pxSize);
        kernel = strel('disk', SeRadius);
        dimensions = size(ChannelContactHeight.Image);
        dimensions2 = size(kernel.Neighborhood);

        % Define kernel center indices
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
                            Cloud_DPs(ii,jj) = ChannelContactHeight.Image(ii,jj)*kernel.Neighborhood(k,l); % Data points within kernel
                            %                         AngleImage(i,j) = AngleImage(i,j) + ChannelHeight.Image(ii,jj)*kernel.Neighborhood(k,l);
                        end
                    end
                end
                %         imagesc(zeros(dimensions(1),dimensions(2))); axis image
                %         hold on
                %         imagesc(Cloud_DPs)
                Cloud_DPs(Cloud_DPs == 0) = NaN;
                [row, col] = ind2sub(size(Cloud_DPs), 1:numel(Cloud_DPs));
                arr = [row(:).*pxSize, col(:).*pxSize, Cloud_DPs(:)]; % Convert axis units to m
                arr(any(isnan(arr), 2), :) = []; % Eliminate rows with nans
                [PlaneAngle] = calculate_angle_topography(arr);
                AngleImage(i,j) = PlaneAngle;
                clear Cloud_DPs
            end
        end
        figure('name', 'Topography angle', 'visible', show_fig);
        imagesc(flipud(AngleImage)); axis image; c = colorbar; 
        c.Location = 'northoutside'; c.Label.String = 'Angle [rad]'; % 90 deg flat 
        set(gca,'FontSize', 16, 'Linewidth', 1.5); axis off
        
        % Angle-based segmentation
        T2 = 1.45; % 83.0788802939694 deg % multithresh(AngleImage); 
        AngleCsBW = imbinarize(AngleImage,T2).*erodedBW2; % Exclude flat areas outside centrosome
        AngleBkgBW = imbinarize(AngleImage,T2).*imcomplement(dilatedBW2); % Glass background area
        figure('name', 'Topography angle segmented','visible',show_fig)
        subplot(1,2,1); hold on
        imagesc(AngleCsBW); axis image; title('centrosome'); axis off 
        set(gca,'FontSize', 16, 'Linewidth', 1.5); 
        subplot(1,2,2); hold on
        imagesc(AngleBkgBW); axis image; title('glass'); axis off 
        set(gca,'FontSize', 16, 'Linewidth', 1.5); 

        % Get indentation modulus Hertz data 
        ChannelEModHertz = E.FM{m}.get_channel(strcat('Indentation Modulus Hertz',s2));
        % From the centrosome's flat region
        CsEModHertz = ChannelEModHertz.Image.*AngleCsBW.*QFits;
        CsEModHertz(CsEModHertz==0)=NaN;
        figure('name', 'Indentation modulus Hertz', 'visible', show_fig); hold on
        imagesc(CsEModHertz*1e-3); axis image; c = colorbar; % kPa
        c.Location = 'northoutside'; c.Label.String = 'Indentation modulus Hertz [kPa]'; 
        axis off; set(gca,'FontSize', 16, 'Linewidth', 1.5); 

        % Calculate centrosome's flat area
        CsFlatProps = regionprops(AngleCsBW, 'Area');
        CsFlatArea = CsFlatProps.Area*pxSize^2;

        % Centrosome height (contact height) and indentation on the flat area
        CsFlatHeight = ChannelContactHeight.Image.*AngleCsBW.*QFits;
        CsFlatHeight(CsFlatHeight==0) = NaN;
        CsFlatMax = max(CsFlatHeight, [], 'all', 'omitnan'); % Maximum height 
        CsFlatPrctile = prctile(CsFlatHeight, 95, "all"); 
        
        ChannelIndenDepth = E.FM{m}.get_channel(strcat('Indentation Depth Hertz',s1,s2));
        CsFlatInden = ChannelIndenDepth.Image.*AngleCsBW.*QFits;
        CsFlatInden(CsFlatInden==0) = NaN;

        % Volume approximation
        CsContactHeight = ChannelContactHeight.Image.*erodedBW2; 
        CsVolume_voxel = CsContactHeight.*pxSize^2; 
        CsVolume = sum(CsVolume_voxel(:)); 
        CsVolumeSphereCap = (pi * CsFlatMax^2 / 3) * (3*CsRadiusXY - CsFlatMax);

        % Save 
        cd(folderPath) 
        filename = strcat('Processed',s2);
        save(filename,'CsArea','CsRadiusXY','CsEModHertz','CsFlatArea','CsFlatHeight','CsFlatMax', 'CsFlatPrctile', 'CsFlatInden','CsVolume', 'CsVolumeSphereCap')

        cd ..
        close all
    end
end