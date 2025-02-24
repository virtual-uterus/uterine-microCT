function WriteSpectrumFile(R,G,B,SpectrumName,FileName)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

RevLabel = {'','reverse'};

fid = fopen(FileName,'w');

fprintf(fid,'gfx create spectrum %s clear overwrite_colour;\n',SpectrumName);

for i=1:size(R,1),
  fprintf(fid,'gfx modify spectrum %s linear %s range %f %f red colour_range %d %d component 1;\n',SpectrumName,RevLabel{R(i,3)+1},R(i,1),R(i,2),min(R(i,3),R(i,4)),max(R(i,3),R(i,4)));
  fprintf(fid,'gfx modify spectrum %s linear %s range %f %f green colour_range %d %d component 1;\n',SpectrumName,RevLabel{G(i,3)+1},G(i,1),G(i,2),min(G(i,3),G(i,4)),max(G(i,3),G(i,4)));
  fprintf(fid,'gfx modify spectrum %s linear %s range %f %f blue colour_range %d %d component 1;\n',SpectrumName,RevLabel{B(i,3)+1},B(i,1),B(i,2),min(B(i,3),B(i,4)),max(B(i,3),B(i,4)));
end;

fclose(fid);

end

