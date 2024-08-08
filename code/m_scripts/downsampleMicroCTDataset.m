function downsampleMicroCTDataset(dir_path, base_name, resize_factor, ...
    varargin)
%DOWNSAMPLEMICROCTDATASET Downsamples a uCT dataset. 
%
%   base_dir is $HOME/Documents/phd/ and set in utils/baseDir()
%
%   Input:
%    - dir_path, path to the directory containing the dataset from base_dir
%    - base_name, name of the dataset.
%    - resize_factor, the factor by which to downsample by, should be < 1.
%    - batch_size, number of image to process in one batch, default 512.
%    - img_extension, extension of the images in the uCT dataset, default
%    bmp.
%   - save_extension, extension of the downsampled images, default png.
%
%   Return:
% Author Alys Clark
% Modified by Emily-Jade Yee
% Modified by Mathias Roesler Nov 2022
narginchk(3, 6);

if nargin < 6
    save_extension = "png";

else
    save_extension = varargin{3};
end

if nargin < 5
    img_extension = "bmp";

else
    img_extension = varargin{2};
end

if nargin < 4
    batch_size = 512;

else
    batch_size = varargin{1};
end

% Directory where images are located
load_directory = join([baseDir(), dir_path, base_name], '/');

save_directory = join([load_directory, "downsampled"], '/');

%% Load and set parameters
% Load properties of the original images that dont change
toml_map = toml.read(join([load_directory, base_name + ".toml"], '/'));
params = toml.map_to_struct(toml_map);
stack_size = double(params.end_nb - params.start_nb);

log_file = join([save_directory, params.prefix + "_downsampled.log"], '/');

if ~isscalar(resize_factor)
    error("The resize factor should be a scalar.")
end

% Estimate new resolution and number of pixels
new_nb_pixel_x = round(params.nb_pixel_x*resize_factor);
new_nb_pixel_y = round(params.nb_pixel_y*resize_factor);
new_stack_size = round(stack_size*resize_factor);
new_resolution_x = (params.nb_pixel_x/new_nb_pixel_x)*params.resolution; % um/pixel x-dir
new_resolution_y = (params.nb_pixel_y/new_nb_pixel_y)*params.resolution; % um/pixel y-dir
new_resolution_z = (stack_size/new_stack_size)*params.resolution; % um/pixel z-dir

nb_runs = ceil(stack_size/batch_size); % Number of times to run loop
img_paths = getImagePaths(load_directory, img_extension);
img_paths = img_paths(params.start_nb:params.end_nb);
img_save_index = 1;

%% Main resizing loop
for run = 1:nb_runs
    disp("Running batch number " + num2str(run) + "/" + num2str(nb_runs));
    first_image_nb = (run-1) * batch_size + 1;
    last_image_nb = run * batch_size;

    if last_image_nb > stack_size
        last_image_nb = stack_size;
    end

    batch_stack_size = last_image_nb - first_image_nb + 1;
    new_batch_stack_size = round(batch_stack_size * resize_factor);

    % Load all images in current batch
    disp("Loading " + num2str(batch_stack_size) + " images in batch");
    img_stack = loadImageStack(img_paths(first_image_nb:last_image_nb));

    % Resize current batch stack
    new_stack = imresize3(uint8(img_stack), resize_factor); 

    % Save current batch stack
    disp("Saving " + num2str(new_batch_stack_size) + " downsampled images");
    saveImageStack(new_stack, save_directory, params.prefix, ...
        img_save_index, save_extension);
    img_save_index = img_save_index + new_batch_stack_size;
end

%% Log required information
file_ID = fopen(log_file, 'w');
fprintf(file_ID, "nb_pixel_x=%3d\n", new_nb_pixel_x);
fprintf(file_ID, "nb_pixel_y=%3d\n", new_nb_pixel_y);
fprintf(file_ID, "stack_size=%3d\n", new_stack_size);
fprintf(file_ID, "pixel_x_res=%03f\n", new_resolution_x);
fprintf(file_ID, "pixel_y_res=%03f\n", new_resolution_y);
fprintf(file_ID, "pixel_z_res=%03f\n", new_resolution_z);
fclose(file_ID);