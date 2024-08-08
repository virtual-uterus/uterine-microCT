function rotated_stack = rotateImageStack(img_stack, region, ...
    nb_used_slices, centreline)
%ROTATEIMAGESTACK Rotates each image in the stack based on the direction of
%the centreline vector
%
%   Input:
%    - img_stack, stack of N MxP images to rotate, img_stack(MxPxN).
%    - region, selects which region to rotate, either left, right.
%    - nb_used_slices, number of slices to take to estimate the centre
%    vector.
%    - centreline, centrepoints for each slice, centreline(6xN).
%
%   Return:
%    - rotated_stack, stack of N MxP rotated images,
%   rotated_stack(MxPxN).
nb_slices = size(img_stack, 3);
rotated_stack = zeros(size(img_stack));

if matches(region, 'left')
    region_nb = 1;
elseif matches(region, 'right')
    region_nb = 5;
else
    error("Error: incorrect value for horn. Should be either left or right")
end

if ~isnumeric(img_stack)
    % Conver to double if not already numeric
    img_stack = double(img_stack);
end

for k = 1:nb_slices

    if mod(k, 50) == 0
        disp("Slice " + k)
    end

    cur_mask = img_stack(:, :, k); % Current mask to rotate

    % Find centre points
    cur_centrepoints = centreline(region_nb:region_nb+1, k)';
    z_centre = nb_used_slices;

    if k < nb_slices-nb_used_slices
        next_centrepoints = centreline(region_nb:region_nb+1, ...
            k + nb_used_slices)';
    else
        % Use previous slices to get rotation vector
        next_centrepoints = centreline(region_nb:region_nb+1, ...
            k - nb_used_slices)';
    end

    % Add the z component
    centre_vector = [next_centrepoints - cur_centrepoints, z_centre];

    % Normalise the centre vector
    centre_vector = centre_vector ./ norm(centre_vector);

    if isequal(centre_vector, [0, 0, 1])
        rotated_stack(:, :, k) = imclose(cur_mask, strel("disk", 1, 4));

    else
        % Create the transformation matrix
        origin = [cur_centrepoints, 0];
        T = findRotationMatrix(centre_vector, [0, 0, 1], origin);

        % Pad the current mask to make a 3D object for rotation
        mask_3D = padarray(cur_mask, [0, 0, 1], 0, 'post');

        % Rotate image and collapse it onto the XY plane
        rotated_mask = sum(imwarp(mask_3D, affine3d(T)), 3);

        % Resize the rotated image to have dimensions of original
        dim_diff = size(cur_mask) - size(rotated_mask);

        if dim_diff(1) > 0
            rotated_mask = padarray(rotated_mask, [dim_diff(1), 0], 'post');

        else
            rotated_mask = rotated_mask(1:end+dim_diff(1), :);
        end

        if dim_diff(2) > 0
            rotated_mask = padarray(rotated_mask, [0, dim_diff(2)], 'post');

        else
            rotated_mask = rotated_mask(:, 1:end+dim_diff(2));
        end

        rotated_stack(:, :, k) = imclose(rotated_mask, strel("disk", 1, 4));

        % Check the ratio of white pixel to remove ill-rotated slices
        nb_w_pixels = sum(sum(rotated_mask > 0));

        if nb_w_pixels / numel(rotated_mask) < 1e-3
            disp("Slice " + k + ": bad rotation");
        end
    end
end