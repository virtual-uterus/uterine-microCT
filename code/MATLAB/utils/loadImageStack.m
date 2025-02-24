function img_stack = loadImageStack(img_paths, mask_paths, xlim, ylim)
%LOADIMAGESTACK Loads the images contains in the img_paths cell array.
%
% If mask_paths is specified, the masks are applied to the images. The
% images are cropped by xlim and ylim.
%   
%   Input:
%    - img_paths, cell array containing the path to every image,
%    img_paths(N x 1).
%    - mask_paths, cell array containing the path to every mask,
%    mask_paths(N x 1), default value [].
%    - xlim, array giving start and end values of x axis crop, default 
%    value [].
%    - ylim, array giving start and end values of y axis crop, default
%    value [].
%
%   Return:
%    - img_stack, matrix containing the images read from the paths. If
%    mask_paths is specified, the masks are applied before return.
if nargin < 4
    ylim = [];
end
if nargin < 3
    xlim = [];
end
if nargin < 2
    mask_paths = [];
end

nb_imgs = length(img_paths);
if nb_imgs ~= length(mask_paths) && ~isempty(mask_paths)
    error("The number of images is different to the number of masks");
end

img = imread(img_paths{1}); % Read first image for size

if isempty(xlim) && isempty(ylim)
    % If no limits provided set x and y limits to be full image
    xlim = [1, size(img, 1)];
    ylim = [1, size(img, 2)];

elseif isempty(xlim)
    % If only x limits not provided set x limits to be full image
    xlim = [1, size(img, 1)];

elseif isempty(ylim)
    % If only y limits not provided set y limits to be full image
    ylim = [1, size(img, 2)];
end

img_stack = zeros(diff(xlim)+1, diff(ylim)+1, nb_imgs, 'uint8');

if ~isempty(mask_paths)
    for k = 1:nb_imgs
        img = imread(img_paths{k}); 
        mask = imread(mask_paths{k});

        if isnumeric(mask)
            % Binarize the image if it is not of logical type
            mask = imbinarize(mask);
        end

        mask = uint8(mask);
        img_stack(:, :, k) = img( ...
            xlim(1):xlim(2), ylim(1):ylim(2)) .* mask( ...
            xlim(1):xlim(2), ylim(1):ylim(2));
    end

else
    for k = 1:nb_imgs
        img = imread(img_paths{k});
        img_stack(:, :, k) = img(xlim(1):xlim(2), ylim(1):ylim(2));
    end

end