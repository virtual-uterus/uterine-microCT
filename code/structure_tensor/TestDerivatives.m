% This script tests different gradient calculations

ImageClass = 'uint8'; MaxVal = 2^8;
Nj = 484; Ni = 162; Nk = 860;
kindex = [1:Nk];
ImageInput = '../Images/SubregionLVFreewall/';
InputFileTemplate = 'SubregionLVFree';
InputFileExtension = 'png';
DigitsInImageSequence = 4; % number of digits in image numbering pattern

% Set the derivative and smoothing template voxel widths
DerivativeTemplateWidth = 5;
%DerivativeTemplateWidth = 3;
SmoothingTemplateWidth = 5;
%SmoothingTemplateWidth = 3;

% Allocate memory
fprintf('... allocating memory ...\n');
I = cast(zeros(Nj,Ni,Nk),ImageClass);

% Load in image set
fprintf('... loading images ...\n');
fstring = sprintf('%s%s%%0%dd.%s',ImageInput,InputFileTemplate,DigitsInImageSequence,InputFileExtension);
for k=1:length(kindex)
  if ~mod(k,100), fprintf(' image: %d\n',k); end
  fnamein = sprintf(fstring,kindex(k)); 
  I(:,:,k) = imread(fnamein);
end

% Pad the image around each edge - reflect image in padding
fprintf('... padding array ...\n');
IPad = padarray(I,[2,2,2],'symmetric');

% convert IPad to a double array between 0 and 1
IPad = double(IPad)/double(max(IPad(:)));

% Use gradient function
[GX,GY,GZ]=gradient(IPad);

% Use other approach with convolution
[Wx,Wy,Wz] = ConstructDerivativeFilters(Nj,Ni,Nk,DerivativeTemplateWidth);
IPadf = fft(reshape(IPad,[prod(size(IPad)),1]));
Wxf = fft(Wx); 
Wyf = fft(Wy); 
Wzf = fft(Wz); 
Di = ifft(IPadf .* Wxf); 
Dj = ifft(IPadf .* Wyf); 
Dk = ifft(IPadf .* Wzf); 

DiR = reshape(Di,size(IPad));
DjR = reshape(Dj,size(IPad));
DkR = reshape(Dk,size(IPad));


% Useful data ranges
PL = 2; % valid derivative padding level
SBi = [1+1*PL:Ni+4-1*PL];
SBj = [1+1*PL:Nj+4-1*PL];
SBk = [1+1*PL:Nk+4-1*PL];

% Raw index locations
[SI0,SJ0,SK0] = ndgrid(1:length(SBi),1:length(SBj),1:length(SBk));

% Level 1
Jii = Di.*Di; Jii = reshape(Jii,[(Ni+4),(Nj+4),(Nk+4)]); 
Jjj = Dj.*Dj; Jjj = reshape(Jjj,[(Ni+4),(Nj+4),(Nk+4)]);
Jkk = Dk.*Dk; Jkk = reshape(Jkk,[(Ni+4),(Nj+4),(Nk+4)]); %clear Dk; 

GX = permute(GX,[2,1,3]);
GY = permute(GY,[2,1,3]);
GZ = permute(GZ,[2,1,3]);

Jii = GY.*GY;
Jjj = GX.*GX;
Jkk = GZ.*GZ;

[Sii1,SI1,SJ1,SK1] = MultigridAveraging(Jii(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
[Sjj1,SI1,SJ1,SK1] = MultigridAveraging(Jjj(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
[Skk1,SI1,SJ1,SK1] = MultigridAveraging(Jkk(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);

% Level 2
[Sii2,SI2,SJ2,SK2] = MultigridAveraging(Sii1,SI1,SJ1,SK1,SmoothingTemplateWidth);
[Sjj2,SI2,SJ2,SK2] = MultigridAveraging(Sjj1,SI1,SJ1,SK1,SmoothingTemplateWidth);
[Skk2,SI2,SJ2,SK2] = MultigridAveraging(Skk1,SI1,SJ1,SK1,SmoothingTemplateWidth);

% Level 3
[Sii3,SI3,SJ3,SK3] = MultigridAveraging(Sii2,SI2,SJ2,SK2,SmoothingTemplateWidth);
[Sjj3,SI3,SJ3,SK3] = MultigridAveraging(Sjj2,SI2,SJ2,SK2,SmoothingTemplateWidth);
[Skk3,SI3,SJ3,SK3] = MultigridAveraging(Skk2,SI2,SJ2,SK2,SmoothingTemplateWidth);

% level 4
[Sii4,SI4,SJ4,SK4] = MultigridAveraging(Sii3,SI3,SJ3,SK3,SmoothingTemplateWidth);
[Sjj4,SI4,SJ4,SK4] = MultigridAveraging(Sjj3,SI3,SJ3,SK3,SmoothingTemplateWidth);
[Skk4,SI4,SJ4,SK4] = MultigridAveraging(Skk3,SI3,SJ3,SK3,SmoothingTemplateWidth);

% level 5
[Sii5,SI5,SJ5,SK5] = MultigridAveraging(Sii4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Sjj5,SI5,SJ5,SK5] = MultigridAveraging(Sjj4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Skk5,SI5,SJ5,SK5] = MultigridAveraging(Skk4,SI4,SJ4,SK4,SmoothingTemplateWidth);

% level 6
[Sii6,SI6,SJ6,SK6] = MultigridAveraging(Sii5,SI5,SJ5,SK5,SmoothingTemplateWidth);
[Sjj6,SI6,SJ6,SK6] = MultigridAveraging(Sjj5,SI5,SJ5,SK5,SmoothingTemplateWidth);
[Skk6,SI6,SJ6,SK6] = MultigridAveraging(Skk5,SI5,SJ5,SK5,SmoothingTemplateWidth);

% level 7
[Sii7,SI7,SJ7,SK7] = MultigridAveraging(Sii6,SI6,SJ6,SK6,SmoothingTemplateWidth);
[Sjj7,SI7,SJ7,SK7] = MultigridAveraging(Sjj6,SI6,SJ6,SK6,SmoothingTemplateWidth);
[Skk7,SI7,SJ7,SK7] = MultigridAveraging(Skk6,SI6,SJ6,SK6,SmoothingTemplateWidth);

% Plots
CI = sub2ind(size(IPad),SJ2(:),SI2(:),SK2(:)); Idx2 = find(IPad(CI) > 1e-6);
CI = sub2ind(size(IPad),SJ3(:),SI3(:),SK3(:)); Idx3 = find(IPad(CI) > 1e-6);
CI = sub2ind(size(IPad),SJ4(:),SI4(:),SK4(:)); Idx4 = find(IPad(CI) > 1e-6);
CI = sub2ind(size(IPad),SJ5(:),SI5(:),SK5(:)); Idx5 = find(IPad(CI) > 1e-6);
CI = sub2ind(size(IPad),SJ6(:),SI6(:),SK6(:)); Idx6 = find(IPad(CI) > 1e-6);
CI = sub2ind(size(IPad),SJ7(:),SI7(:),SK7(:)); Idx7 = find(IPad(CI) > 1e-6);
figure(4); clf; plot([2,3,4,5,6,7],[median(Skk2(Idx2)),median(Skk3(Idx3)),median(Skk4(Idx4)),median(Skk5(Idx5)),median(Skk6(Idx6)),median(Skk7(Idx7))],'r-*',[2,3,4,5,6,7],[median(Sjj2(Idx2)),median(Sjj3(Idx3)),median(Sjj4(Idx4)),median(Sjj5(Idx5)),median(Sjj6(Idx6)),median(Sjj7(Idx7))],'g-*',[2,3,4,5,6,7],[median(Sii2(Idx2)),median(Sii3(Idx3)),median(Sii4(Idx4)),median(Sii5(Idx5)),median(Sii6(Idx6)),median(Sii7(Idx7))],'b-*');
