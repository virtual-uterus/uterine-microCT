function angle = ComputeFibreAngle(fibre, centrepoints, cur_X, cur_Z, ...
    nb_used_slices)
%COMPUTEFIBREANGLE Computes the angle in degrees of the fibre relative to
%the plane given by the vector between two centre points. 
%
%   Inputs:
%    - fibre, 3D array containing the XYZ coordinates of the fibre.
%    - centrepoints, 6xN, list of centrepoints.
%    - cur_X, x coordinate of the current fibre in the general coordinate
%    system.
%    - cur_Z, z coordinate of the current fibre in the general coordinate
%    system.
%    - nb_used_slices, number of slices to use to estimate centre vector.
%
%   Return:
%    - angle, angle between the fibre and the current plane in degrees.
if isempty(centrepoints)
    z_vector = [0; 0; 1];

else
    cur_Z = round(cur_Z);
    cur_X = round(cur_X);
    
    if all(centrepoints(3:4, cur_Z))
        % If a middle point is found
        % Define the z vector as the centre vector between current slice
        % and the current + 5th slice.
        cur_idx = 3;
    else
        % If there are left and right centre points find the nearest to
        % cur_X
        [~, cur_idx] = min(abs(centrepoints(1:4:end, cur_Z) - cur_X));

        if cur_idx == 2
            % Get the starting index of the last centre point
            cur_idx = 5;
        end
    end

    cur_centrepoint = centrepoints(cur_idx:cur_idx+1, cur_Z);

    if cur_Z + nb_used_slices > size(centrepoints, 2)
        % If there are not enough slices use the previous slices
        % instead
        next_centrepoint = centrepoints(cur_idx:cur_idx+1, cur_Z - nb_used_slices);
        z_vector = cur_centrepoint - next_centrepoint;
        z_vector = [z_vector; -nb_used_slices]; % Append the z component

    else
        next_centrepoint = centrepoints(cur_idx:cur_idx+1, cur_Z + nb_used_slices);
        z_vector = next_centrepoint - cur_centrepoint;
        z_vector = [z_vector; nb_used_slices]; % Append the z component

    end
end

angle = rad2deg(acos(fibre * (z_vector/norm(z_vector))));
end