function WriteOrthoFile(file_name, fibres, sheets, normals, fibre_angles)
%WRITEORTHOFILE Writes the ortho file based on fibres, sheets, and normals
%vectors.
%
%   Args:
%    - file_name, path to the ortho file to write.
%    - fibres, fibre orientation vectors fibres(N, 3).
%    - sheets, sheet orientation vectors sheets(N, 3).
%    - normals, normal orientation vectors normals(N, 3).
%    - fibre_angles, angles relative to the centreline fibre_angles(N, 1).
%
%   Return:
if nargin < 5
    fibre_angles = [];
end
%% Check the fibres, sheets, and normals have the same size
if ~all(size(fibres) == size(sheets))
    error("Error: fibres and sheets should have the same size.")
end
if ~all(size(fibres) == size(normals))
    error("Error: fibres and normals should have the same size.")
end

fid = fopen(file_name,'w');
fprintf(fid,'%d\n', length(fibres));  % Print number of elements

for i = 1:length(fibres)
    fprintf(fid, '%f %f %f', fibres(i, 1), fibres(i, 2), fibres(i, 3));
    fprintf(fid, ' %f %f %f ', sheets(i, 1), sheets(i, 2), sheets(i, 3));
    fprintf(fid, '%f %f %f ', normals(i, 1), normals(i, 2), normals(i, 3));
    if length(fibre_angles) ~= 0
        fprintf(fid, '%f\n', fibre_angles(i));
    end
end

fclose(fid);
end

