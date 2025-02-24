% This script finds the principal axes of the chamber masks and stores them
% for subsequent calculations

% Heart 1 images
MaskPathFull = '../ImagesHeart1/ChamberMask/';
OutputPathReduced = '../DataHeart1/STBinary/';
MaskPrefix = 'ChamberMask_';
LVGroup = 2; RVGroup = 1;
Njm = 2048; Nim = 2352; Nkm = 1493;
MRange = [1:1493];

% Read in mask data and 
kstart = MRange(1); kend = MRange(Nkm); 
I3D = false(Njm,Nim,Nkm); 
for k=kstart:kend
  fnamein = sprintf('%s%s%04d.png',MaskPathFull,MaskPrefix,k);
  M = ~(imread(fnamein) == 0);
  I3D(:,:,k-kstart+1) = M;
end
I3D = permute(I3D,[2,1,3]);

% Find two largest volumes
fprintf(' ... finding largest volumes ...\n');
%I3D = logical(I3D);
stats = regionprops(I3D,'Area','PixelList','Centroid');
[Order,Oloc]=sort(cat(1,stats(:).Area),'descend');

% 3. Use PCA to find main axis
fprintf('... find main axis with PCA \n');
CoefsLV = pca(stats(LVGroup).PixelList);
COMLV = stats(LVGroup).Centroid;
%COMLV = COMLV([2,1,3]);

CoefsRV = pca(stats(RVGroup).PixelList);
COMRV = stats(RVGroup).Centroid;
%COMRV = COMRV([2,1,3]);

% 4.1. Write data to mat file for helix angle computation
save(sprintf('%sPermuteChamberAxisDataForHelixAngleCalcs.mat',OutputPath),'CoefsLV','CoefsRV','COMLV','COMRV');

 
