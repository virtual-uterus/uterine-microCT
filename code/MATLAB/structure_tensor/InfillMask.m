%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% InfillMask.m
%
% This script uses key frames to fill in  mask removal regions
% The filled in removal mask is dilated in each plane and its complement
% applied to the original mask.
%
% Written by: Mark Trew, November 2021.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: Load in key frames and morph masks between them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Images from Fetal Heart 1
% Nj = 2048; Ni = 2352;
% kindex = [1521:1751];
% fullkindex = [260:1752];
% KeyStep = 10;
% ImageInput = '../ImagesHeart1/MaskSampleEdit/';
% FileTemplate = 'FH1_PTA_20_7_21_';
% ImageOutput = '../ImagesHeart1/MaskSampleFill/';
% MaskInput = '../ImagesHeart1/MaskCombined/';
% MaskedInput = '../ImagesHeart1/MaskedCombined/';
% MaskOutput = '../ImagesHeart1/MaskEdited/';
% MaskedOutput = '../ImagesHeart1/MaskedEdited/';

% Images from Air dry H1C1H
% Nj = 1589; Ni = 1771;
% kindex = [1401:1689];
% fullkindex = [370:1689];
% KeyStep = 32;
% ImageInput = '../Images_H1C1H/KeyImages/';
% FileTemplate = 'H1C1H_';
% ImageOutput = '../Images_H1C1H/KeyImagesFilled/';
% MaskInput = '../Images_H1C1H/Mask/';
% MaskedInput = '../Images_H1C1H/Masked/';
% MaskOutput = '../Images_H1C1H/MaskClean/';
% MaskedOutput = '../Images_H1C1H/MaskedClean/';

% Images from Air dry H1C1H
Nj = 1613; Ni = 2118;
kindex = [1416:1616];
fullkindex = [1:1630];
KeyStep = 10;
ImageInput = '../Images_H1C1H/KeyImages/';
FileTemplate = 'H1C1H_';
ImageOutput = '../Images_H1C1H/KeyImagesFilled/';
MaskInput = '../Images_H1C1H/Mask/';
MaskedInput = '../Images_H1C1H/Masked/';
MaskOutput = '../Images_H1C1H/MaskClean/';
MaskedOutput = '../Images_H1C1H/MaskedClean/';

InputFileExtension = 'png';
DigitsInImageSequence = 4; % number of digits in image numbering pattern

% Allocate memory
KeyFrameIndex = [kindex(1):KeyStep:kindex(end)];
KeyFrames = false(Nj,Ni,length(KeyFrameIndex)); 
InfillFrames = false(Nj,Ni,length(kindex));

% Read key frames
fprintf(' ... Reading key frames and morphing ...\n');
for k=1:length(KeyFrameIndex)
    Temp = imread(sprintf('%s%s%04d.png',ImageInput,FileTemplate,KeyFrameIndex(k)));
    KeyFrames(:,:,k) = Temp(:,:,1);
end

% Morph between key frames
k = 1;
for kk=1:length(KeyFrameIndex)-1
    for kifill=1:KeyStep
        InfillFrames(:,:,k) = MaskMorph(KeyFrames(:,:,kk),KeyFrames(:,:,kk+1),KeyFrameIndex(kk),KeyFrameIndex(kk+1),kindex(k));
        k = k+1;
    end
end
InfillFrames(:,:,end) = KeyFrames(:,:,end);

% Write morphed infilled images
fprintf(' ... Writing morphed images ...\n');
for k=1:length(kindex)
    imwrite(InfillFrames(:,:,k),sprintf('%s%s%04d.png',ImageOutput,FileTemplate,kindex(k)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2: Apply morphed masks to full masks and masked images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in mask and masked tissue
Mask = false(Nj,Ni,length(fullkindex));
Tissue = uint8(false(Nj,Ni,length(fullkindex)));
% MK = [1150:1661];

[dummy,IntersectIndex,dummy]=intersect(fullkindex,kindex);
[dummy,DifferenceIndex]=setdiff(fullkindex,kindex);

fprintf(' ... reading in images ...\n');
for k=1:length(fullkindex)
    Mask(:,:,k) = imread(sprintf('%s%s%04d.png',MaskInput,FileTemplate,fullkindex(k)));
    Tissue(:,:,k) = imread(sprintf('%s%s%04d.png',MaskedInput,FileTemplate,fullkindex(k)));
end

% Apply masking
fprintf(' ... applying masking ...\n');
for k=1:length(IntersectIndex)
    Mask(:,:,IntersectIndex(k)) = ~imdilate(InfillFrames(:,:,k),ones(3,3)).*Mask(:,:,IntersectIndex(k));
end

% Only keep largest volume
fprintf(' ... finding largest volume ...\n');
stats = regionprops(Mask,'Area','PixelIdxList');
[mv,mloc]=max(cat(1,stats(:).Area));
Mask = false(Nj,Ni,length(fullkindex));
Mask(stats(mloc).PixelIdxList) = 1;

% Write mask files
fprintf(' ... Writing mask and masked tissue ...\n');
for k=1:length(fullkindex)
    imwrite(Mask(:,:,k),sprintf('%s%s%04d.png',MaskOutput,FileTemplate,fullkindex(k)),'BitDepth',1);
    imwrite(uint8(Mask(:,:,k)).*Tissue(:,:,k),sprintf('%s%s%04d.png',MaskedOutput,FileTemplate,fullkindex(k)));
end
