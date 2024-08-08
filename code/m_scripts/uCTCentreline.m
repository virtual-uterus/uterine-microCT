function uCTCentreline(dir_path, base_name, regions, downsampled, ST, ...
    extension)
%UCTCENTRELINE Computes the centreline for the dataset provided by
%base_name given the selected regions. 
%
%   If the region is "both", the centreline is smoothed using a
%   Savitzky-Golay method.
%   
%   base_dir is $HOME/Documents/phd/ and set in utils/baseDir()
%
%   Input:
%    - dir_path, path to the directory containing the dataset from base_dir
%    - base_name, name of the dataset.
%    - region, either left, right or both, used to sort the centrepoints. 
%    - downsampled, true if the dataset has been downsampled, default value
%    is true.
%    - ST, true if the dataset is located in the ST folder, otherwise
%    located in the muscle_segmentation folder, default value is true.
%    - extension, extension of the images to load, default value is png.
%
%   Return: 
if nargin < 6
    extension = "png";
end
if nargin < 5
    ST = true;
end
if nargin < 4
    downsampled = true;
end

% Directory where images are located
load_directory = join([baseDir(), dir_path, base_name], '/');

if downsampled
    % If using the downsampled dataset
    load_directory = join([load_directory, "downsampled"], '/');
end

if ST
    % Deal with final location
    load_directory = join([load_directory, "ST/mask"], '/');
else
    load_directory = join([load_directory, "muscle_segmentation"], '/');
end

for k = 1:length(regions)
    region = regions(k);
    disp("Processing region: " + region)

    tmp_load_directory = load_directory;

    if strcmp(region, "both")
        mask_paths = getImagePaths(tmp_load_directory, extension);
    else
        tmp_load_directory = join([load_directory, region], '/');
        mask_paths = getImagePaths(tmp_load_directory, extension);
    end

    mask_stack = loadImageStack(mask_paths);

    nb_slices = size(mask_stack, 3);
    centreline = zeros(6, nb_slices); % Placeholder for 3 centre points

    for l = 1:nb_slices
        centre_points = findCentrepoints(mask_stack(:, :, l), region);
        centreline(:, l) = reshape(centre_points', [6, 1]);
    end

    if strcmp(region, "both")
        for m = 1:size(centreline, 1)
            % Smooth centreline coordinates
            end_idx = find(centreline(m, :), 1, 'last');
            centreline(m, 1:end_idx) = smoothdata(centreline(m, 1:end_idx), ...
                2, "sgolay");
        end
    end

    disp("Saving centreline")
    save(tmp_load_directory + "/centreline.mat", "centreline");

    clear centreline
end
