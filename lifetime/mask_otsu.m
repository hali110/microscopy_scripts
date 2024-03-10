function [mask_fill, centroid] = mask_otsu(int_shape, thresh, disk)

%     figure
%     rawIm=int_shape;
%     rawImNorm = rawIm/nanmax(rawIm(:));
%     
%     [OT,EM] = graythresh(rawImNorm); % OT: Otsu Threshold, EM: Effectiveness Metric, 0 = 1 population, 1 = 2 populations
%     
%     BWs = imbinarize(rawImNorm); % Otsu
%     
%     se = strel('disk',3); % r = 3 works well with Otsu thresholding
%     BWsdil = imdilate(BWs,se);
%     BWdfill = imfill(BWsdil,'holes');
% 
%     BWdfill = imfill(BWs, 'holes');
% 
%     % Wnobord = imclearborder(BWdfill,4);
%     BWfinal = BWdfill;
%     
%     mask=NaN(512);
%     
%     % mask(BWs)=1;
%     mask(BWfinal)=1;
%     % mask(~BWfinal) = nan;
%     centroid=1;

    n_thresholds = thresh;
    disk_size = disk;

    level = multithresh(int_shape,n_thresholds);
    mask = imquantize(int_shape,level);
    % figure, imagesc(mask),
    mask_bw = mask > 1;
    % figure, imagesc(mask_bw)
    
    SE = strel('disk',disk_size);
    mask_bw1 = imclose(mask_bw,SE);
    % figure, imagesc(mask_bw1)
    mask_fill = imfill(mask_bw1, 'holes');

    % Calculate properties of all connected components
    L = bwlabeln(mask_fill);
    props = regionprops(L, 'Area', 'Perimeter', 'Eccentricity', 'Centroid');
    areas = cat(1,props.Area);
    ind = find(areas == max(areas));

    % figure(1), imagesc(mask_fill)
    % figure(2), imagesc(int_shape)

    centroid = props(ind).Centroid;

    mask_fill = double(mask_fill);
    mask_fill(~mask_fill) = NaN;

end