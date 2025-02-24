function processed_img_stack = preprocessImageStack(img_stack, fct_nb, ...
    varargin)
%PREPROCESSIMAGE Preprocesses the inputed image based on the selected
%method.
%   Inputs:
%    - img, image to preprocess.
%    - fct_nb, number associated with the preprocessing function:
%     1 -> imadjust
%     2 -> histeq
%     3 -> imsharpen
%    - varargin, extra arguments for the imflatfield function.
narginchk(2, 4);

if nargin == 3
    img_strel = strel('disk', varargin{1});
end

processed_img_stack = uint8(zeros(size(img_stack)));

parfor k = 1:size(processed_img_stack, 3)
    switch fct_nb
        case 1
            processed_img_stack(:, :, k) = imadjust(img_stack(:, :, k));
            processed_img_stack(:, :, k) = imclose( ...
                processed_img_stack(:, :, k), img_strel);
        case 2
            processed_img_stack(:, :, k) = histeq(img_stack(:, :, k));
            processed_img_stack(:, :, k) = imclose( ...
                processed_img_stack(:, :, k), img_strel);
        case 3
            processed_img_stack(:, :, k) = imsharpen(img_stack(:, :, k));
            processed_img_stack(:, :, k) = imclose( ...
                processed_img_stack(:, :, k), img_strel);    
    end

end