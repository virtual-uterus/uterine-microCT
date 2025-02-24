function mask_stack = shapeSegmentation(img_stack, threshold, ...
    morph_size, neighborhood_size)
%TISSUESEGMENTATION Segments the shape from the background of the images
%in the provided stack.
%
%   Input:
%    - img_stack, stack of images to segment.
%    - threshold, threshold value for the binary segmentation.
%    - morph_size, size of the morphological operation structuring element.
%    - neighborhood_size, size of the neighborhood of the median filter.
%
%   Return:
%    - mask_stack, stack of the masks associated with each image.
mask_stack = zeros(size(img_stack));
morph_SE = ones(morph_size, morph_size);
medfilt_neighborhood = [neighborhood_size, neighborhood_size];

parfor k = 1:size(img_stack, 3)
    img = img_stack(:, :, k);
    mask = zeros(size(img));

    mask(img > threshold) = 255;
    mask = imopen(imclose(mask, morph_SE), morph_SE);
    mask = imfill(mask, 'holes');
    mask_stack(:, :, k) = medfilt2(mask, medfilt_neighborhood);

end
end