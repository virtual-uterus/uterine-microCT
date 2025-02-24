function removeOvaries(dir_path, base_name, downsampled, extension, start_nb)
%REMOVEOVARIES Removes the ovaries from the segmentation masks of a uCT
%dataset.
%
%   base_dir is $HOME/Documents/phd/ and set in utils/baseDir()
%
%   Input:
%    - dir_path, path to the directory containing the dataset from base_dir
%    - base_name, name of the dataset.
%    - downsampled, true if the dataset has been downsampled, default value
%    is true.
%    - extension, extension of the segmentation masks, default png.
%    - start_nb, number at which to start saving images, default 1.
%   Return:
if nargin < 6
    start_nb = 1;
end
if nargin < 5
    extension = "png";
end
if nargin < 4
    downsampled = true;
end

% Directory where images are located
load_directory = join([baseDir(), dir_path, base_name], '/');

if downsampled
    % If using the downsampled dataset
    load_directory = join([load_directory, "downsampled"], '/');
    toml_map = toml.read(join([load_directory, base_name + "_downsampled.toml"], '/'));
else
    % Use the non-downsampled TOML file
    toml_map = toml.read(join([load_directory, base_name + ".toml"], '/'));
end

% Load parameters
params = toml.map_to_struct(toml_map);

img_paths = getImagePaths(load_directory, extension);

disp('Loading image stack')
img_stack = loadImageStack(img_paths);

middle_pixel = floor(max(size(img_stack)) / 2);
[~, largest_dim] = max([size(img_stack, 1), size(img_stack, 2)]);
nb_img = size(img_stack, 3);

disp('Removing left ovary')
for k = params.left_ovary:nb_img
    % Remove ovary on the left side of the image
    switch largest_dim
        case 1
            img_stack(:, 1:middle_pixel, k) = 0;

        case 2
            img_stack(1:middle_pixel, :, k) = 0;
    end
end

disp('Removing right ovary')
for k = params.right_ovary:nb_img
    % Remove ovary on the right side of the image
    switch largest_dim
        case 1
            img_stack(:, middle_pixel:end, k) = 0;

        case 2
            img_stack(middle_pixel:end, :, k) = 0;
    end
end

disp('Saving image stack')
saveImageStack(img_stack, load_directory, params.prefix, start_nb, extension);

end