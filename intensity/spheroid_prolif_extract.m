close all, clear, clc 

%% Define parameters to read .nd2 files

% Define source and target folders and add bioformats to path

origin_folder = pwd;
target_folder = '/Volumes/TMR/tmr_data/Nikon A1R/2021-10-21';

%% Read .nd2 files and convert them to .mat files

% Get filenames
cd(target_folder)
files = dir('*.nd2');
filenames = cell(size(files,1),1);
for i = 1:size(files,1)
    filenames(i) = cellstr(files(i).name); 
end





for i = 1:length(filenames)

    % Make directory to save min intensity projections
    % mkdir(files(i).name(end-22:end-21))

    % Access reader without loading .nd2 file
    reader = bfGetReader(filenames{i});

        
    % Extract metadata
    omeMeta = reader.getMetadataStore();
    
    % X-axis data
    Pixel_n.X = omeMeta.getPixelsSizeX(0).getValue();                                                   % No. of pixels in x-direction
    Pixel_um.X = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER).doubleValue();     % Size of each pixel in um in x-direction
    % Y-axis data
    Pixel_n.Y = omeMeta.getPixelsSizeY(0).getValue();                                                   % No. of pixels in y-direction
    Pixel_um.Y = omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER).doubleValue();     % Size of each pixel in um in y-direction
    % Z-axis data
    Pixel_n.Z = omeMeta.getPixelsSizeZ(0).getValue();                                                   % No. of pixels in z-direction
    Pixel_um.Z = omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER).doubleValue();     % Size of each pixel in um in z-direction
    % Time data
    Pixel_n.T = omeMeta.getPixelsSizeT(0).getValue();
    
    
    data = zeros(reader.getSeriesCount() - 1, Pixel_n.Y,Pixel_n.X,Pixel_n.Z);
    for series = 0:reader.getSeriesCount() - 1 
        reader.setSeries(series)

        % Loop through each Time, Z, and Channel and save all planes
        for z = 0:Pixel_n.Z - 1

            % Find what plane the specified z,c,t are, (t,c = 0,0 for yap - t,c = 0,1 for taz,cntrl)
            plane = reader.getIndex(z, 1, 0) + 1;

            % Read said plane and save
            data(series+1,:,:,z+1) = bfGetPlane(reader, plane);
                    
        end
        
        disp(['saved series:' num2str(series)])


       
        
    end
    % Take a minimum intenisty and save it
    DIC_tmp = min(data, [], 4);
    save([files(i).folder '/analysis/DIC_tmp.mat'], 'DIC_tmp')

        
    
end
    
    

%% Binary Mask
clc, close all, clearvars -except filenames target_folder origin_folder DIC_tmp

% Initialize parameters and variables
crop = 210;
n_thresholds = 4;
disk_size = 4;
reader = bfGetReader(filenames{1});

files = {'DIC_tmp.mat'};

% Loop through and threshold
    

for i = 1:length(files)


    reader = bfGetReader(filenames{i});
    omeMeta = reader.getMetadataStore();
    um_per_pix = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER).doubleValue();
    spheroid_rad_um = NaN(reader.getSeriesCount(), 1);
    % Load data from mat file
    disp(['Loading:' files{i}])
    % load(['analysis/' files{i}])

    for series = 1:size(DIC_tmp, 1)

        frame(:,:) = DIC_tmp(series, crop:end-crop, crop:end-crop);
        % frame(:,:) = DIC_tmp(series,:,:);
        % figure, imagesc(frame)
        
        fig = figure('Position', [150 300 1100 500]);
        subplot(1, 2, 1), imagesc(frame)
    
        % Binary mask using Otsu's method
        level = multithresh(frame,n_thresholds);
        mask = imquantize(frame,level);
        % figure, imagesc(mask),
        mask_bw = mask < n_thresholds;
        % mask_bw = mask == 1;
        % figure, imagesc(mask_bw)
        
        SE = strel('disk',disk_size);
        mask_bw1 = imclose(mask_bw,SE);
        % figure, imagesc(mask_bw1)

        imq = imfill(mask_bw1, 'holes');
        % figure, imagesc(imq)
        
        
        % Calculate properties of all connected components
        L = bwlabeln(imq);
        props = regionprops(L, 'Area', 'Perimeter', 'Eccentricity', 'Centroid');
        
        % Extract the radius and centroid of the element having the largest area
        areas = cat(1,props.Area);
        ind = find(areas == max(areas));
        spheroid_rad_pix = (props(ind).Area / pi) ^ (1/2);
        spheroid_rad_um(series) =  spheroid_rad_pix * um_per_pix;

        subplot(1, 2, 2), imagesc(imq)
        hold on
        plot(props(ind).Centroid(1), props(ind).Centroid(2), 'ro', 'MarkerSize', 30);

        
        % clear level mask mask_bw SE mask_bw1 L props areas ind
    end

    % radii2save{i} = spheroid_rad_um;


end
    


% save([target_folder '/analysis/radii.mat'], 'spheroid_rad_um')
    
    
