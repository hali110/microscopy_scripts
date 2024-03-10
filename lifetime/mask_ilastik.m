function [mask_fill, centroid] = mask_ilastik(mask_fold,imag_name, ilastik_params)
%%
% Make masks for FLIM data: take the Ilastik multi-layer segmentation 
% ouput and create a binary mask.
% Mask code: 
% c = Cytoplasm
% n = Nuclei
% b = Background
%%

c = ilastik_params(1);
n = ilastik_params(2);
b = ilastik_params(3);



% initialize images
mask = [];
Cyto_Mask = [];
Nucl_Mask = []; 
Back_Mask = [];

% import ilastik map generated via pixel classification
mask = imread([mask_fold '/' imag_name '-Intensity_Simple Segmentation.tif']);

% collapse the masks to projections
Cyto_Mask = any(mask == c, 3);
Nucl_Mask = any(mask == n, 3);
Back_Mask = any(mask == b, 3);

% Calculate properties of all connected components
L = bwlabeln(Cyto_Mask);
props = regionprops(L, 'Area', 'Perimeter', 'Eccentricity', 'Centroid');
areas = cat(1,props.Area);
ind = find(areas == max(areas));
if isempty(ind)
    centroid = NaN;
else
    centroid = props(ind).Centroid;
end

mask_fill = double(Cyto_Mask);
mask_fill(~Cyto_Mask) = NaN;

end