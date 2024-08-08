function T = findRotationMatrix(u, v, origin)
%FINDROTATIONMATRIX Finds the rotation matrix between two vectors.
%
%   The rotation axis is the cross product of u and v.
%   Input:
%    - u, first vector u(3x1).
%    - v, second vector v(3x1).
%    - origin, origin of the translation.
%
%   Return:
%    - T, transformation matrix in homogeneous coordinates T(4x4)
cross_product = cross(u, v);

% Normalize the rotation axis
cross_product = cross_product / norm(cross_product);

% Calculate the dot product between u and v
dot_product = dot(u, v);

% Calculate the skew-symmetric matrix
K = [0, -cross_product(3), cross_product(2);
     cross_product(3), 0, -cross_product(1);
     -cross_product(2), cross_product(1), 0];

% Calculate the rotation matrix
R = eye(3) + K + K^2 * (1 - dot_product) / norm(cross_product)^2;

% Create the transformation matrix in homogeneous coordinates
T = eye(4);
T(1:3, 1:3) = R;
T(4, 1:3) = origin;
end