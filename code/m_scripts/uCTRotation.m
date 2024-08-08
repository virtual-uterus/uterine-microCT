function uCTRotation(dir_path, base_name, regions, type, downsampled, ...
    extension)
%UCTROTATION Computes the centreline for the dataset provided by
%base_name given the selected regions.
%   
%   base_dir is $HOME/Documents/phd/ and set in utils/baseDir()
%
%   Input:
%    - dir_path, path to the directory containing the dataset from base_dir
%    - base_name, name of the dataset.
%    - region, either left, right or ["left", "right"], used to sort the centrepoints. 
%    - type, segmentation type, {fat, tissue, shape, muscle}, default value
%    is muscle.
%    - downsampled, true if the dataset has been downsampled, default value
%    is true.
%    - extension, extension of the images to load, default value is png.
%
%   Return: 
if nargin < 6
    extension = "png";
end
if nargin < 5
    downsampled = true;
end
if nargin < 4
    type = "muscle";
end

% Directory where images are located
load_directory = join([baseDir(), dir_path, base_name], '/');

if downsampled
    % If using the downsampled dataset
    load_directory = join([load_directory, "downsampled"], '/');
    toml_map = toml.read(join([load_directory, ...
        base_name + "_downsampled.toml"], '/'));
else
    % Use the non-downsampled TOML file
    toml_map = toml.read(join([load_directory, base_name + ".toml"], '/'));
end

% Add the segmentation type 
load_directory = join([load_directory, type + "_segmentation/"], '/');

% Load parameters
params = toml.map_to_struct(toml_map);
start_nb = params.thickness.start_nb;
nb_used_slices = double(params.nb_used_slices);

mask_paths = getImagePaths(load_directory, extension);
mask_stack = loadImageStack(mask_paths);
load(load_directory + 'centreline.mat', 'centreline');

for k = 1:length(regions)
    region = regions(k);

    if strcmp(region, "left")
        end_nb = params.thickness.left.end_nb;
    elseif strcmp(region, "right")
        end_nb = params.thickness.right.end_nb;
    else
        error("Error: invalid horn selection.");
    end

    disp("Rotating region: " + region);
    rotated_stack = rotateImageStack( ...
        mask_stack(:, :, start_nb:end_nb), region, nb_used_slices, ...
        centreline(:, start_nb:end_nb)); 

    disp("Saving region: " + region);
    save_directory = join([load_directory, region], '/');
    saveImageStack(rotated_stack, save_directory, ...
        params.prefix, start_nb, extension);

    clear rotated_stack % Save memory
end
