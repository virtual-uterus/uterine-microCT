function WriteOrthoFile(file_name, fibres, sheets, normals)
%WRITEORTHOFILE Writes the ortho file based on fibres, sheets, and normals
%vectors.
%
%   Args:
%    - file_name, path to the ortho file to write.
%    - fibres, fibre orientation vectors fibres(N, 3).
%    - sheets, sheet orientation vectors fibres(N, 3).
%    - normals, normal orientation vectors fibres(N, 3).
%
%   Return:
%% Check the fibres, sheets, and normals have the same size
if ~all(size(fibres) == size(sheets))
    error("Error: fibres and sheets should have the same size.")
end
if ~all(size(fibres) == size(normals))
    error("Error: fibres and normals should have the same size.")
end

fid = fopen(file_name,'w');
fprintf(fid,'%d', length(fibres));  % Print number of elements

for i = 1:length(fibres)
    fprintf(fid, '%f %f %f', fibres(i, 1), fibres(i, 2), fibres(i, 3));
    fprintf(fid, ' %f %f %f', sheets(i, 1), sheets(i, 2), sheets(i, 3));
    fprintf(fid, ' %f %f %f', normals(i, 1), normals(i, 2), normals(i, 3));
end

fclose(fid);
end

