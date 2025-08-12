function mask_stack = endometriumSegmentation(img_stack, ...
    muscle_segmentation_stack)
%ENDOMETRIUMSEGMENTATION Segments the endometrium from the background of 
%the images in the provided stack.
%
%   Input:
%    - img_stack, stack of images to segment.
%    - muscle_segmentation_stack, stack of mask of the muscle segementation
%
%   Return:
%    - mask_stack, stack of the masks associated with each image.
mask_stack = zeros(size(img_stack));


for k = 1:size(img_stack, 3)
    img = img_stack(:, :, k);
    muscle_mask = muscle_segmentation_stack(:, :, k);
    filled_mask = imfill(muscle_mask, 'holes');
    inner_mask = filled_mask - muscle_mask;
    mask_stack(:, :, k) = imbinarize(img .* (inner_mask / 255));
end
end