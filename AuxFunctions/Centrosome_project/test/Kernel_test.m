format long g;
format compact;
workspace;  % make sure the workspace panel is showing

close all
clc
% clear

% image1 = imread('cameraman.tif'); 
% image1 = double(image1);
image1 = ChanHeight.Image*1e9; 
% image1 = image1';
kernel = strel('disk', 5);

% function [filtered] = basic_convolution(image,kernel)
dimensions = size(image1);
dimensions2 = size(kernel.Neighborhood);

% define kernel center indices
kernelCenter_x = round(dimensions2(1)/2);
kernelCenter_y = round(dimensions2(2)/2);

image2 = zeros(dimensions(1),dimensions(2));
for i = 1:dimensions(1)
    for j = 1:dimensions(2)
        for k = 1:dimensions2(1)
            for l = 1:dimensions2(2)
                ii = i+(k-kernelCenter_x);
                jj = j+(l-kernelCenter_y);
                if (ii >= 1 && ii <= dimensions(1) && jj >= 1 && jj <= dimensions(2))
                    Cloud_DPs(ii,jj) = image1(ii,jj)*kernel.Neighborhood(k,l); % data points within kernel
                    %                         image2(i,j) = image2(i,j) + image(ii,jj)*kernel.Neighborhood(k,l);
                end
            end
        end
%         imagesc(zeros(dimensions(1),dimensions(2))); axis image
%         hold on
%         imagesc(Cloud_DPs)
        Cloud_DPs(Cloud_DPs == 0) = NaN; 
        [row, col] = ind2sub(size(Cloud_DPs), 1:numel(Cloud_DPs));
        arr = [row(:), col(:), Cloud_DPs(:)];
        arr(any(isnan(arr), 2), :) = [];% eliminate rows with nans
        [PlaneAngle] = calculate_topography_angle(arr); 
        image2(i,j) = PlaneAngle;
        clear Cloud_DPs
    end
%     filtered = image2;
end

imagesc(image2); axis image; colorbar
figure()
imagesc(image1); axis image; colorbar

figure()
test_filt = stdfilt(image1, kernel.Neighborhood);
imagesc(test_filt); axis image; colorbar

segmented_test = image2; 
% segmented_test(segmented_test>90.1) = 0; 
% segmented_test(segmented_test<89.9) = 0; 
figure()
hold on 
imagesc(segmented_test); axis image 