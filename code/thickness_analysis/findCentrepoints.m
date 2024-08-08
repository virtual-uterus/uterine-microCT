function centrepoints = findCentrepoints(mask, region)
%FINDCENTREPOINTS Finds the centrepoint of a mask given the region.
%
%   If there is no clear separation between the left and right horn, three
%   points are given: left, centre, and right.
%
%   Input:
%    - mask, binary mask.
%    - region, place the centrepoint in the corresponding region if only
%    one centrepoint is found.
%
%   Return:
%    - centrepoints, coordinates of the centrepoints, the first pair of
%    is the left point, the second pair the middle point, the third pair
%    the right point, centrepoints(3, 2).
if ~islogical(mask)
    % Convert to logical if not because regionprops only works on logical
    mask = imbinarize(mask);
end

proper_size = 300; % Define how big a region is to be considered
centrepoints = zeros(3, 2); % Placeholder for the centrepoints

filled_mask = imfill(mask, "holes");
centre_regions = and(not(mask), filled_mask);

filled_props = regionprops(filled_mask, 'Area', 'Centroid');
proper_region = sum([filled_props.Area] > proper_size);

centre_props = regionprops(centre_regions, 'Area', 'Centroid');
nb_holes = sum([centre_props.Area] > proper_size);

if proper_region == 1 && nb_holes <= 1
    % Use the number of holes to filter out the region. If there are two
    % then in the body, if there is 1 then single horn

    % Single horn and can exit early
    centrepoint = centre_props( ...
        [centre_props.Area] > proper_size).Centroid;

    if strcmp(region, "left")
        % Left horn
        centrepoints(1, :) = centrepoint;

    else
        % Right horn
        centrepoints(3, :) = centrepoint;
    end

    return;
end

% Create arrays to recuperate region properties if not exited early
areas = zeros(length(centre_props), 1);
centroids = zeros(length(centre_props), 2);

for k = 1:length(centre_props)
    areas(k) = centre_props(k).Area;
    centroids(k, :) = centre_props(k).Centroid;
end

[~, ind] = maxk(areas, 2); % Find two largest areas indices
centroids = centroids(ind, :); % Get the correct centroids

% Refine the centres by using the skeleton
skeleton = bwskel(centre_regions);
[idx_y, idx_x] = find(skeleton == 1);

for k = 1:size(centroids, 1)
    differences = [idx_x, idx_y] - centroids(k, :);
    [~, min_idx] = min(vecnorm(differences'));
    centroids(k, :) = [idx_x(min_idx), idx_y(min_idx)];
end

if size(centroids, 1) == 1
    % If only on centroid sort it based on region
    if strcmp(region, "left")
        % Left horn
        centrepoints(1, :) = centroids;

    else
        % Right horn
        centrepoints(3, :) = centroids;
    end

else
    % Find the left and right centroids and sort them
    if centroids(1, 1) > centroids(2, 1)
        centrepoints(1, :) = centroids(2, :);
        centrepoints(3, :) = centroids(1, :);

    else
        centrepoints(1, :) = centroids(1, :);
        centrepoints(3, :) = centroids(2, :);
    end
end

if proper_region == 1 && nb_holes > 1
    % Use the number of holes to filter out the region. If there are two
    % then in the body, if there is 1 then single horn
    % In the body to find middle centre point

    u = centrepoints(1, :) - centrepoints(3, :); % Orientation vector
    u = u ./ norm(u); % Normalised vector
    
    t = -max(size(mask)):max(size(mask));
    x_points = centrepoints(1, 1) + u(1) .* t;
    y_points = centrepoints(1, 2) + u(2) .* t;

    % Filter out points that are not in the image
    valid_indices = (x_points >= min(centrepoints(1, 1), centrepoints(3, 1))) ...
        & (x_points <= max(centrepoints(1, 1), centrepoints(3, 1))) & ...
        (y_points >= min(centrepoints(1, 2), centrepoints(3, 2))) & ... 
        (y_points <= max(centrepoints(1, 2), centrepoints(3, 2)));
    x_points = x_points(valid_indices);
    y_points = y_points(valid_indices);

    % Convert subscripts to linear indices to find central region
    mask_idx = sub2ind(size(mask), round(y_points), round(x_points));
    centre_region = find(mask(mask_idx) == 1);

    % Get the centre of the region and add it to the centrepoints
    middle_centrepoint = mask_idx(centre_region( ...
        round(length(centre_region) / 2)));
    [centrepoints(2, 2), centrepoints(2, 1)] = ind2sub( ...
        size(mask), middle_centrepoint);

end