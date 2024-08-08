function mask_stack = fatSegmentation(img_stack, morph_size)
%FATSEGMENTATION Segments the fat from the uterine horn from each image in
%the provided stack.
%
%   Input:
%    - img_stack, stack of images to segment.
%    - morph_size, size of the morphological operation structuring element.
%
%   Return:
%    - mask_stack, stack of the masks associated with each image.
mask_stack = zeros(size(img_stack));

% Structuring elements
morph_SE = strel('disk', morph_size);

%% Image segmentation with k-means
parfor k = 1:size(img_stack, 3)
    mask = mask_stack(:, :, k);

    [idx, centers] = imsegkmeans(img_stack(:, :, k), 3);
    [~, max_ind] = max(centers); % The desired layer is associated with max

    layers_idx = idx == max_ind;
    mask(layers_idx) = 1;
    mask = imopen(mask, morph_SE);
    mask_stack(:, :, k) = imbinarize(mask);
end
end