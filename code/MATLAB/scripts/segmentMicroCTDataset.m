function segmentMicroCTDataset(dir_path, base_name, segmentation_type, ...
    downsampled, extension, start_nb)
%SEGMENTMICROCTDATASET Segments a uCT dataset based on the segmentation 
%type. 
%
%   base_dir is $HOME/Documents/phd/ and set in utils/baseDir()
%
%   Input:
%    - dir_path, path to the directory containing the dataset from base_dir
%    - base_name, name of the dataset.
%    - segmentation_type, type of segmentation to use:
%      1 is tissue segmentation
%      2 is muscle layer segmentation
%      3 is fat segmentation
%      4 is shape segmentation
%    - downsampled, true if the dataset has been downsampled, default value
%    is true.
%    - extension, extension of the segmentation masks, default png.
%    - start_nb, number at which to start saving images, default 1.
%   Return:
if nargin < 7
    start_nb = 1;
end
if nargin < 6
    extension = "png";
end
if nargin < 5
    downsampled = true;
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

% Load parameters
params = toml.map_to_struct(toml_map);
preprocess = params.preprocess;
morph_size = double(params.morph_size);

% Get image paths
img_paths = getImagePaths(load_directory, extension);

disp('Loading image stack')
img_stack = loadImageStack(img_paths);

disp('Spliting image stack')
stack_cell = splitImageStack(img_stack, params.split_nb);
mask_cell = cell(size(stack_cell)); % Empty cell array for the masks

for k = 1:length(stack_cell)
    img_stack = stack_cell{k};
    disp('Processing stack ' + string(k))
    if preprocess ~= 0
        disp('    Preprocessing image stack')
        img_stack = preprocessImageStack(img_stack, preprocess, ...
            morph_size);
    end

    disp('    Segmenting image stack')
    switch segmentation_type
        case 1
            % Tissue segmentation
            save_dir= join([load_directory, "tissue_segmentation"], '/');
            mask_stack = tissueSegmentation(img_stack, ...
                params.segmentation.tissue.threshold, ...
                params.segmentation.tissue.morph_size, ...
                params.segmentation.tissue.neighborhood_size);

        case 2
            % Muscle segmentation
            save_dir = join([load_directory, "muscle_segmentation"], '/');
            mask_stack = muscleSegmentation(img_stack, ...
                params.segmentation.muscle.morph_size, ...
                params.segmentation.muscle.neighborhood_size);

        case 3
            % Fat segmentation
            save_dir= join([load_directory, "fat_segmentation"], '/');
            mask_stack = fatSegmentation(img_stack, ...
                params.segmentation.fat.morph_size);

        case 4
            % Shape segmentation
            save_dir= join([load_directory, "shape_segmentation"], '/');
            mask_stack = shapeSegmentation(img_stack, ...
                params.segmentation.shape.threshold, ...
                params.segmentation.shape.morph_size, ...
                params.segmentation.shape.neighborhood_size);
    end

    mask_cell{k} = mask_stack;
end

disp('Fusing mask stack')
mask_stack = fuseImageStacks(mask_cell);

disp('Saving mask stack')
saveImageStack(mask_stack, save_dir, params.prefix, start_nb, extension);


end