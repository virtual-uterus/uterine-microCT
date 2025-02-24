%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script determines a ventricle COM and long axis vector required for 
% computing helix angles across the ventricles.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Load in mask images

% Lower ventricles Heart 2
% MaskPath = '../Images/Mask/';
% MaskPrefix = 'FH2_PTA_12_8_21_';
% OutputPath = '../Data/STBinary/';
% Njm = 2048; Nim = 2048; Nkm = 1024;
% MRange = 141+[100:1123];

% Upper ventricles Heart 2
% MaskPath = '../Images/Mask_Part2_CleanMask/AllFramesMasked/';
% MaskPrefix = 'FH2_PTA_12_8_21_';
% Njm = 2048; Nim = 2048; Nkm = 512;
% MRange = [1150:1661];

% Full Heart 1
% MaskPath = '../ImagesHeart1/MaskedEdited/';
% MaskPrefix = 'FH1_PTA_20_7_21_';
% OutputPath = '../DataHeart1/STBinary/';
% ImageOutputPath = '../ImagesHeart1/MaskCutHelixSurf_8bit_512x512x256/';
% Njm = 2048; Nim = 2352; Nkm = 1493;
% MRange = [260:1752];

% Air Dry Heart 2
% MaskPath = '../ImagesHeartDry/Mask/';
% MaskPrefix = 'AirDry_';
% OutputPath = '../DataHeartDry/STBinary/';
% ImageOutputPath = '../ImagesHeartDry/HelixAngleImages/DisplayMask/';
% Njm = 1237; Nim = 1237;
% MRange = [548:1060];
% Nkm = length(MRange);

% Air Dry H1C1H
% MaskPath = '../Images_H1C1H/Mask/';
% MaskPrefix = 'H1C1H_';
% OutputPath = '../Data_H1C1H/STBinary/';
% ImageOutputPath = '../ImagesHeartDry/HelixAngleImages/DisplayMask/';
% Njm = 1589; Nim = 1771;
% MRange = [370:1689];
% Nkm = length(MRange);

% Air Dry H1FGR1H
% MaskPath = '../Images_H1FGR1H/Mask/';
% MaskPrefix = 'H1FGR1H_';
% OutputPath = '../Images_H1FGR1H/Data_H1FGR1H/STBinary/';
% ImageOutputPath = '../Images_H1FGR1H/HelixAngleImages/DisplayMask/';
% Njm = 1586; Nim = 2083;
% MRange = [1:1460];
% Nkm = length(MRange);

% Air Dry H1C1H 20 Sept 2022
MaskPath = '../Images_H1C1H/Mask/';
MaskPrefix = 'H1C1H_';
OutputPath = '../Images_H1C1H/Data_H1C1H/STBinary/';
ImageOutputPath = '../Images_H1C1H/HelixAngleImages/DisplayMask/';
Njm = 1613; Nim = 2118;
MRange = [1:1616];
Nkm = length(MRange);

% Read in mask data
kstart = MRange(1); kend = MRange(Nkm); 
I3D = false(Njm,Nim,Nkm); 
for k=kstart:kend
  fnamein = sprintf('%s%s%04d.png',MaskPath,MaskPrefix,k);
  M = ~(imread(fnamein) == 0);
  I3D(:,:,k-kstart+1) = M;
end

% 2. Find mask voxel locations - note, should this be the ventricular wall
% or the LV lumen?
fprintf('... find voxel locations \n');
stats = regionprops(I3D,'Area','PixelList','Centroid');
[mv,mloc]=max(cat(1,stats(:).Area));
PixelLocs = stats(mloc).PixelList;

% 3. Use PCA to find main axis of voxel cloud
fprintf('... find main axis with PCA \n');
Coefs = pca(PixelLocs);
PA = Coefs(:,3);
PA = PA([2,1,3]);
R = zeros(3,3);
RComp = [0,-1;1,0];
[mv,midx]=min(abs(PA));
idx = setdiff([1:3],midx);
R(idx,idx) = RComp;
r = R*PA;
r = r/norm(r);
q = cross(PA,r);

% 4. Find the centre of mass
COM = stats(mloc).Centroid;
COM = COM([2,1,3]);

% 4.1. Write data to mat file for helix angle computation
save(sprintf('%sAxisDataForHelixAngleCalcs.mat',OutputPath),'Coefs','COM');

return;

% 5. Export visualisation data - subsampled point cloud (data points), COM (data point), long axis (1D elements),
% finite element plane orthogonal to long axis (face in 3D space).
% 5.1. 
DIdx = 1:20000:size(PixelLocs,1);
WriteGeneralExdataFile(PixelLocs(DIdx,1),PixelLocs(DIdx,2),PixelLocs(DIdx,3),[1:length(DIdx)]',[],sprintf('%sMaskPoints',OutputPath),'MaskPoints',{});
% sculp mask
DIdx = 1:1000:size(PixelLocs,1);
DIdx = DIdx(find(PixelLocs(DIdx,3) < round(COM(3)) | (PixelLocs(DIdx,1) > round(COM(1)) & PixelLocs(DIdx,2) > round(COM(2)))));
WriteGeneralExdataFile(PixelLocs(DIdx,1),PixelLocs(DIdx,2),PixelLocs(DIdx,3),200000+[1:length(DIdx)]',[],sprintf('%sMaskPointsSculpt',OutputPath),'MaskPointsSculpt',{});

% 5.2. 
XYZA = [COM-700*PA';COM;COM+500*PA'];
WriteGeneralLineExnodeExelemFile(XYZA(:,1),XYZA(:,2),XYZA(:,3),[1,2;2,3],[11:13]',[11:12]',[],sprintf('%sHeartAxis',OutputPath),'HeartAxis',{});

% 5.3
D = 700;
XYZP = [COM-D*r'-D*q';COM+D*r'-D*q';COM-D*r'+D*q';COM+D*r'+D*q'];
WriteGeneralPlaneExnodeExelemFile(XYZP(:,1),XYZP(:,2),XYZP(:,3),[1,2,3,4],[14:17]',[13],[],sprintf('%sHeartHelixPlane',OutputPath),'HeartHelixPlane',{});

% 5.4 Cut lower part of heart below the COM
IdxBelow = find(PA(1)*PixelLocs(DIdx,1)+PA(2)*PixelLocs(DIdx,2)+PA(3)*PixelLocs(DIdx,3) < (dot(PA,COM)-10));
WriteGeneralExdataFile(PixelLocs(DIdx(IdxBelow),1),PixelLocs(DIdx(IdxBelow),2),PixelLocs(DIdx(IdxBelow),3),[1:length(DIdx(IdxBelow))]',[],sprintf('%sMaskPointsBelow',OutputPath),'MaskPoints',{});

IdxBelow = find(PA(2)*PixelLocs(:,1)+PA(1)*PixelLocs(:,2)+PA(3)*PixelLocs(:,3) < (dot(PA,COM)-10));
CutMask = false(size(I3D));
PtIdx = sub2ind(size(I3D),PixelLocs(IdxBelow,2),PixelLocs(IdxBelow,1),PixelLocs(IdxBelow,3));
CutMask(PtIdx) = 1;
CutMaskR = imresize3(CutMask,[512,512,256],'Antialiasing',true);
fprintf('... writing images ...\n');
for k=1:256
     fnameout = sprintf('%sMask_%04d.png',ImageOutputPath,k);
     imwrite(CutMaskR(:,:,k),fnameout);
%     imwrite(imrotate(flipud(CutMaskR(:,:,k)),-90),fnameout);
end

