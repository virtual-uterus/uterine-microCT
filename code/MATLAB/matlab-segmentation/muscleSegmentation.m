function mask_stack = muscleSegmentation(img_stack, morph_size, ...
    neighborhood_size)
%MUSCLESEGMENTATION Segments the muscle layers from each image in the
%provided stack.
%
%   Input:
%    - img_stack, stack of images to segment.
%    - morph_size, size of the morphological operation structuring element.
%    - neighborhood_size, size of the neighborhood of the median filter.
%
%   Return:
%    - mask_stack, stack of the masks associated with each image.
mask_stack = zeros(size(img_stack));
morph_SE = ones(morph_size, morph_size);
medfilt_neighborhood = [neighborhood_size, neighborhood_size];

%% Image segmentation with k-means
parfor k = 1:size(img_stack, 3)
    img = img_stack(:, :, k);
    mask = zeros(size(img));

    if ~isequal(img, mask)
        [idx, centers] = imsegkmeans(img, 3);
        [~, max_ind] = max(centers); % The desired layer is associated with the max

        layers_idx = idx == max_ind;

        mask(layers_idx) = 1;

        mask = medfilt2(mask, medfilt_neighborhood);
        mask = imopen(imclose(mask, morph_SE), morph_SE);
    end

    mask_stack(:, :, k) = mask;
end

end