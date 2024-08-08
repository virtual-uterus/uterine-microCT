%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MaskChambers.m
%
% This script loads in mask files and segments out the ventricles.
%
% It downsizes the masks for display purposes, finds the left ventricle. 
%
% Helix plane is determined from the LV lumen principal axis.
%
% Updated by: Mark Trew, June 2022.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Data and path locations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Myocardium images - lower ventricles
% InputPath = '../Data/STBinary/';
% OutputPath = '../Images/ChamberMask_2048x2048x1024/';
% MaskPath = '../Images/Mask/';
% MaskPrefix = 'FH2_PTA_12_8_21_';
% OutputPrefix = 'ChamberMask_';
% Njm = 2048; Nim = 2048; Nkm = 1024;
% MRange = 141+[100:1123];

% Myocardium images - upper ventricles
% OutputPath = '../Images/Mask_Part2_CleanMask/ChamberMask_2048x2048x512/';
% MaskPath = '../Images/Mask_Part2_CleanMask/AllFramesMasked/';
% MaskPrefix = 'FH2_PTA_12_8_21_';
% OutputPrefix = 'ChamberMask_';
% Njm = 2048; Nim = 2048; Nkm = 512;
% MRange = [1150:1661];
 
% Heart 1 images
% MaskPath = '../ImagesHeart1/MaskEdited/';
% OutputPathFull = '../ImagesHeart1/ChamberMask/';
% OutputPathReduced = '../ImagesHeart1/ResizedDisplay/ChamberMask_8bit_512x588x373/';
% MaskPrefix = 'FH1_PTA_20_7_21_';
% OutputPrefix = 'ChamberMask_';
% Njm = 2048; Nim = 2352; Nkm = 1493;
% MRange = [260:1752];
% NewSize = [512,588,373];

% Air Dry H1C1H
% MaskPath = '../Images_H1C1H/MaskClean/';
% OutputMaskPathReduced = '../Images_H1C1H/MaskClean_8bit_512x571x425/';
% OutputPathFull = '../Images_H1C1H/ChamberMask/';
% OutputPathReduced = '../Images_H1C1H/ChamberMask_8bit_512x571x425/';
% DataOutputPath = '../Data_H1C1H/STBinary/';
% MaskPrefix = 'H1C1H_';
% OutputPrefix = 'H1C1H_ChamberMask_';
% Njm = 1589; Nim = 1771;
% MRange = [370:1689];
% Nkm = length(MRange);
% NewSize = [512,571,425];

% Air Dry H1FGR1H
% MaskPath = '../Images_H1FGR1H/Mask/';
% OutputMaskPathReduced = '../Images_H1FGR1H/Mask_8bit_512x521x396/';
% OutputPathFull = '../Images_H1FGR1H/ChamberMask/';
% OutputPathReduced = '../Images_H1FGR1H/ChamberMask_8bit_512x521x396/';
% DataOutputPath = '../Images_H1FGR1H/Data_H1FGR1H/STBinary/';
% MaskPrefix = 'H1FGR1H_';
% OutputPrefix = 'H1FGR1H_ChamberMask_';
% Njm = 1586; Nim = 2083;
% MRange = [1:1460];
% Nkm = length(MRange);
% NewSize = [512,521,396];

% Air Dry H1C1H 20 Sept 2022
MaskPath = '../Images_H1C1H/Mask/';
OutputMaskPathReduced = '../Images_H1C1H/Mask_8bit_512x530x403/';
OutputPathFull = '../Images_H1C1H/ChamberMask/';
OutputPathReduced = '../Images_H1C1H/ChamberMask_8bit_512x530x403/';
DataOutputPath = '../Images_H1C1H/Data_H1C1H/STBinary/';
MaskPrefix = 'H1C1H_';
OutputPrefix = 'H1C1H_ChamberMask_';
Njm = 1613; Nim = 2118;
MRange = [1:1616];
Nkm = length(MRange);
NewSize = [512,530,403];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load mask data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' ... Loading image data ... \n');

I3D = false(Njm,Nim,Nkm); 
parfor k=1:Nkm
  fnamein = sprintf('%s%s%04d.png',MaskPath,MaskPrefix,MRange(k));
  M = ~(imread(fnamein) == 0);
  I3D(:,:,k) = M;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fill ventricles.
%
% Mask needs to have a zero upper surface so that isosurfaces are correctly
% calculated. - only on 8 bit
%
% For H1C1H need to use a two step approach.
% 1. 2D hole fill of top surface
% 2. 3D hole fill of ventricles
% 3. Subtract original mask from filled mask - this will find the LV.
% 4. Find largest volume - this will be the LV
% 5. Do imclose with ones(9,9,9)
% 6. 3D hole fill of ventricles
% 7. Subtract first filled mask from second filled mask - this will find
% the RV
% 8. Find largest volume - this will be the RV
%
% Otherwise
% 1. 2D hole fill of top surface
% 2. 3D fill of ventricles
% 3. Subtract orginal mask from filled mask
% 4. Find two largest volumes (LV largest, RV second largest)
%
% Alternative - not exact but should work.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('... Iterative 2D filling ...\n');
Working = I3D;
% Working(:,:,end) = imfill(Working(:,:,end),'holes'); %1
% Working1 = imfill(Working,'holes'); %2
% Working2 = imsubtract(Working1,I3D); %3
% stats = regionprops(Working2,'Area','PixelIdxList'); %4
% [Order,Oloc]=sort(cat(1,stats(:).Area),'descend'); %4
% LVMask = false(Njm,Nim,Nkm); %4
% LVMask(stats(Oloc(1)).PixelIdxList) = 1; %4
% %RVMask(stats(Oloc(2)).PixelIdxList) = 1; %4 - not for H1C1H
% 
% % Additional for H1C1H
% Working3 = imclose(Working1,ones(9,9,9)); %5
% Working3 = imfill(Working3,'holes'); %6
% Working3 = imsubtract(Working3,Working1); %7

% Alternative - iterate and do 2D fill
parfor k=1:Nkm
    Working(:,:,k) = logical(imfill(Working(:,:,k),'holes'));
end
Working = logical(imsubtract(Working,I3D));

fprintf('... Cleaning ventricular mask ...\n');
stats = regionprops(Working,'Area','PixelIdxList');
[Order,Oloc]=sort(cat(1,stats(:).Area),'descend');
LVMask = false(Njm,Nim,Nkm);
LVMask(stats(Oloc(1)).PixelIdxList) = 1; 
RVMask = false(Njm,Nim,Nkm);
RVMask(stats(Oloc(2)).PixelIdxList) = 1; 
Mask = LVMask+RVMask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Write full chamber mask and resized masks
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write full chamber mask to file
fprintf('...Writing full chamber mask...\n');
parfor k=1:size(Mask,3)
     fnameout = sprintf('%s%s%04d.png',OutputPathFull,OutputPrefix,k);
     imwrite(Mask(:,:,k),fnameout);
end

% Resize chamber mask and convert to 8 bit
fprintf('...Resize chamber mask...\n');
MaskR = imresize3(255*uint8(Mask),NewSize,'Antialiasing',true);
MaskR(:,:,end) = 0;

% Write resized chamber mask
fprintf('...Writing resized chamber mask...\n');
parfor k=1:size(MaskR,3)
     fnameout = sprintf('%s%s%04d.png',OutputPathReduced,OutputPrefix,k);
     imwrite(MaskR(:,:,k),fnameout);
end

% Resize ventricle tissue mask and convert to 8 bit
fprintf('...Resize ventricular tissue mask...\n');
VMaskR = imresize3(255*uint8(I3D),NewSize,'Antialiasing',true);
VMaskR(:,:,end) = 0;

% Write resized mask
fprintf('...Writing resized ventricular tissue mask...\n');
parfor k=1:size(VMaskR,3)
     fnameout = sprintf('%s%s%04d.png',OutputMaskPathReduced,OutputPrefix,k);
     imwrite(VMaskR(:,:,k),fnameout);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute ventricle chamber measures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Left Ventricle
LVSurf = imsubtract(LVMask,imerode(LVMask,ones(3,3,3)));
LVSurf(:,:,end-1:end) = 0;
% Right Ventricle
RVSurf = imsubtract(RVMask,imerode(RVMask,ones(3,3,3)));
RVSurf(:,:,end-1:end) = 0;

% cumulated surface to volume ratios
PartialSum_LV_Surf = zeros(Nkm,1);
PartialSum_LV_Vol = zeros(Nkm,1);
PartialSum_RV_Surf = zeros(Nkm,1);
PartialSum_RV_Vol = zeros(Nkm,1);

parfor k=1:Nkm
    Slice = LVSurf(:,:,k);
    PartialSum_LV_Surf(k) = sum(Slice(:));
    Slice = LVMask(:,:,k);
    PartialSum_LV_Vol(k) = sum(Slice(:));
    Slice = RVSurf(:,:,k);
    PartialSum_RV_Surf(k) = sum(Slice(:));
    Slice = RVMask(:,:,k);
    PartialSum_RV_Vol(k) = sum(Slice(:));
end

CSSurfLV = cumsum(PartialSum_LV_Surf);
CSVolLV = cumsum(PartialSum_LV_Vol);
CSSurfRV = cumsum(PartialSum_RV_Surf);
CSVolRV = cumsum(PartialSum_RV_Vol);

% Save data in files
save(sprintf('%sCumulativeSVandEAData.mat',DataOutputPath),'CSSurfLV','CSVolLV','CSSurfRV','CSVolRV','PartialSum_LV_Surf','PartialSum_LV_Vol','PartialSum_RV_Vol');

% Plots - for visual assessmetn of surface to volume (cumulative 3D) and
% edge length to cross section area (2D)
% %figure(2); h=plot([1:Nkm]/1000,CSSurfLV./(CSVolLV+1e06),'r',[1:Nkm]/1000,CSSurfRV./(CSVolRV+1e06),'b','LineWidth',1.5);
% figure(2); h=plot([1:length(CSSurfLV)]/1000,CSSurfLV(:)./(CSVolLV(:)+1e06),'r',[1:length(CSSurfRV)]/1000,CSSurfRV(:)./(CSVolRV(:)+1e06),'b','LineWidth',1.5); axis([0 1.5 0 0.1]);
% legend(h,'Left ventricle','Right ventricle','FontName','Arial','FontSize',12);
% ylabel('Surface to volume ratio (1/mm)','FontName','Arial','FontSize',12);
% set(gca,'FontName','Arial','FontSize',12);
% xlabel('Cumulative distance (mm)','FontName','Arial','FontSize',12);
% %print('../DataHeart1/STBinary/CumulatedSurfaceVolumeRatio.png','-dpng','-r1200')
% print('../Data/STBinary/CumulatedSurfaceVolumeRatio.png','-dpng','-r1200')
% 
% %figure(3); h=plot([1:Nkm-2]/1000,PartialSum_LV_Surf(1:end-2)./(PartialSum_LV_Vol(1:end-2)+1e06),'r',[1:Nkm-2]/1000,PartialSum_RV_Surf(1:end-2)./(PartialSum_RV_Vol(1:end-2)+1e06),'b','LineWidth',1.5);
% figure(3); h=plot([1:length(PartialSum_LV_Surf)]/1000,PartialSum_LV_Surf(:)./(PartialSum_LV_Vol(:)+1e06),'r',[1:length(PartialSum_RV_Surf)]/1000,PartialSum_RV_Surf(:)./(PartialSum_RV_Vol(:)+1e06),'b','LineWidth',1.5);axis([0 1.5 0 0.025]);
% legend(h,'Left ventricle','Right ventricle','FontName','Arial','FontSize',12);
% xlabel('Distance (mm)','FontName','Arial','FontSize',12);
% ylabel('Edge length to slice area ratio (1/mm)','FontName','Arial','FontSize',12);
% set(gca,'FontName','Arial','FontSize',12);
% %print('../DataHeart1/STBinary/SampledEdgeLengthVolumeRatio.png','-dpng','-r1200')
% print('../Data/STBinary/SampledEdgeLengthVolumeRatio.png','-dpng','-r1200')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute LV principal component for helix angle.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find mask voxel locations - note, should this be the ventricular wall
% or the LV lumen?
fprintf('... Find LV lumen voxel locations ... \n');
stats = regionprops(LVMask,'Area','PixelList','Centroid');
[mv,mloc]=max(cat(1,stats(:).Area));
PixelLocs = stats(mloc).PixelList;

% Use PCA to find main axis of voxel cloud
fprintf('... find main axis with PCA \n');
Coefs = pca(PixelLocs);
% PA = Coefs(:,3);
% PA = PA([2,1,3]);
PA = Coefs(:,1);
R = zeros(3,3);
RComp = [0,-1;1,0];
[mv,midx]=min(abs(PA));
idx = setdiff([1:3],midx);
R(idx,idx) = RComp;
r = R*PA;
r = r/norm(r);
q = cross(PA,r);

% Find the centre of mass
COM = stats(mloc).Centroid;
%COM = COM([2,1,3]);

% Write data to mat file for helix angle computation
save(sprintf('%sLVAxisDataForHelixAngleCalcs.mat',DataOutputPath),'Coefs','COM');

% Write visualisation files
% Write axis file
D = 500;
XYZA = [COM-D*PA';COM;COM+D*PA'];
%WriteGeneralLineExnodeExelemFile(XYZA(:,1),XYZA(:,2),XYZA(:,3),[1,2;2,3],[11:13]',[11:12]',[],sprintf('%sHeartAxis',DataOutputPath),'HeartAxis',{});
WriteGeneralLineExnodeExelemFile(XYZA(:,2),XYZA(:,1),XYZA(:,3),[1,2;2,3],[11:13]',[11:12]',[],sprintf('%sHeartAxis',DataOutputPath),'HeartAxis',{});

% Write plane file at the COM 
D = 700;
XYZP = [COM-D*r'-D*q';COM+D*r'-D*q';COM-D*r'+D*q';COM+D*r'+D*q'];
WriteGeneralPlaneExnodeExelemFile(XYZP(:,1),XYZP(:,2),XYZP(:,3),[1,2,3,4],[14:17]',[13],[],sprintf('%sHeartHelixPlane',DataOutputPath),'HeartHelixPlane',{});
