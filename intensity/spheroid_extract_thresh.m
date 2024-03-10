close all, clear, clc 

%% Define parameters to read .nd2 files

% Define source and target folders and add bioformats to path

origin_folder = pwd;
target_folder = '/Volumes/NikonAX/Data/2022-04-19';

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
    mkdir(files(i).name(end-7:end-4))

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
    
    
    
    for series = 0:reader.getSeriesCount() - 1 
        data = zeros(Pixel_n.Y,Pixel_n.X,Pixel_n.Z,Pixel_n.T,reader.getSizeC());
        reader.setSeries(series)

        % Loop through each Time, Z, and Channel and save all planes
        for t = 0:Pixel_n.T - 1
            for z = 0:Pixel_n.Z - 1
                for c = 0:reader.getSizeC() - 1

                    % Find what plane the specified z,c,t are
                    plane = reader.getIndex(z, c, t) + 1;

                    % Read said plane and save
                    data(:,:,z+1,t+1,c+1) = bfGetPlane(reader, plane);
                    disp(['saved ztc:' num2str(z) num2str(t) num2str(c)])
                    
                end
            end
        end
        


        % Take a minimum intenisty and save it
        DIC_tmp = min(data(:,:,:,:,2), [], 3);
        save([files(i).folder '/' files(i).name(end-7:end-4) '/series' num2str(series) '.mat'], 'DIC_tmp')

        
        
    end
    
end
    
    

%% Binary Mask

% Initialize parameters and variables
crop = 80;
n_thresholds = 1;
disk_size = 8;
reader = bfGetReader(filenames{i});
omeMeta = reader.getMetadataStore();
um_per_pix = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER).doubleValue();

spheroid_rad_um = NaN(length(filenames), reader.getSeriesCount());


% Loop through and threshold
for day = 1:length(filenames)
    
    % Find all .mat files
    series_data = dir(['day' files(day).name(end-4) '/*.mat']);

    for sample = 1:length(series_data)
        
        % Load data from mat file
        disp(['Loading:' series_data(sample).folder '/' series_data(sample).name])
        load([series_data(sample).folder '/' series_data(sample).name])
        % frame(:,:) = DIC_tmp(crop:end-crop, crop:end-crop, 1, end);
        frame = DIC_tmp(:,:,1,end);
        figure, imagesc(frame)

        % Binary mask using Otsu's method
        level = multithresh(frame,n_thresholds);
        mask = imquantize(frame,level);
        % figure, imagesc(mask),
        mask_bw = imcomplement(mask > n_thresholds);
        % figure, imagesc(mask_bw)
        
        SE = strel('disk',disk_size);
        mask_bw1 = imclose(mask_bw,SE);
        figure, imagesc(mask_bw1)
        
        % Calculate properties of all connected components
        L = bwlabeln(mask_bw1);
        props = regionprops(L, 'Area', 'Perimeter', 'Eccentricity', 'Centroid');
        
        % Extract the radius and centroid of the element having the largest area
        areas = cat(1,props.Area);
        ind = find(areas == max(areas));
        spheroid_rad_pix = (props(ind).Area / pi) ^ (1/2);
        spheroid_rad_um(day, sample) =  spheroid_rad_pix * um_per_pix;
        
        % clear level mask mask_bw SE mask_bw1 L props areas ind
    end
    
end

% save([target_folder '/radii.mat'], 'spheroid_rad_um')
    
    
