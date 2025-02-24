function [img_stack, mask_stack] = readNiftiSegmentation(dir_path, ...
    base_name, type, downsampled, binarise)
%READNIFTISEGMENTATION Reads the segmentation and images from the nifti
%files. The masks are applied to the images.
%   
%   base_dir is $HOME/Documents/phd/ and set in utils/baseDir()
%
%   Input:
%    - dir_path, path to the directory containing the dataset from base_dir
%    - base_name, name of the dataset.
%    - type, segmentation type, {fat, tissue, shape, muscle}, default value
%    is muscle.
%    - downsampled, true if the dataset has been downsampled, default value
%    is true.
%    - binarise, true if the masks should be binarise, defautl value is
%    true. 
%
%   Return:
%    - img_stack, uint8, stack of images read from the nifti file. The
%    segmentation masks have been applied.
%    - mask_stack, logical or uint8, stack of segmentation masks.
if nargin < 3
    type = "muscle";
end
if nargin < 4
    downsampled = true;
end
if nargin < 5
    binarise = true;
end

% Directory where images are located
load_directory = join([baseDir(), dir_path, base_name], '/');

if downsampled
    % If using the downsampled dataset
    load_directory = join([load_directory, "downsampled"], '/');
end

mask_file = base_name + '_' + type + "_segmentation.nii.gz";
img_file = base_name + ".nii.gz";

mask_stack = niftiread(join([load_directory, mask_file], '/'));
img_stack = niftiread(join([load_directory, img_file], '/'));

% Binarize the masks
if binarise
    mask_stack = imbinarize(mask_stack); 

else
    mask_stack = uint8(mask_stack);
end

% Permute the first two columns of images and masks
mask_stack = permute(mask_stack, [2 1 3]);
img_stack = permute(img_stack, [2, 1, 3]);

% Apply masks
img_stack = uint8(mask_stack) .* img_stack;

end