% This script resizes the ventricular mask - and erodes it slightly - and
% the masked images so they can be used for display purposes.

% Heart 1 images
InputMaskPath = '../ImagesHeart1/MaskEdited/';
InputTissuePath = '../ImagesHeart1/MaskedEdited/';
OutputMaskPath = '../ImagesHeart1/ResizedDisplay/MaskErode_8bit_512x588x373/';
OutputMaskPath2 = '../ImagesHeart1/ResizedDisplay/MaskErode_Cut_8bit_512x588x373/';
OutputMaskPath3 = '../ImagesHeart1/ResizedDisplay/MaskErode_Cut_Opposite_8bit_512x588x373/';
OutputTissuePath = '../ImagesHeart1/ResizedDisplay/Masked_Cut_8bit_512x588x373/';
InputMaskPrefix = 'FH1_PTA_20_7_21_';
InputTissuePrefix = 'FH1_PTA_20_7_21_';
OutputMaskPrefix = 'Mask_';
OutputTissuePrefix = 'Tissue_';
Njm = 2048; Nim = 2352; Nkm = 1493;
MRange = [260:1752];
NewSize = [512,588,373];

% Heart 2 images
% InputMaskPath = '../Images/Mask/';
% InputTissuePath = '../Images/Masked/';
% OutputMaskPath = '../Images/MaskErode_8bit_512x512x256/';
% OutputMaskPath2 = '../Images/MaskErode_Cut_8bit_512x512x256/';
% OutputTissuePath = '../Images/Masked_8bit_512x512x256/';
% InputMaskPrefix = 'FH2_PTA_12_8_21_';
% InputTissuePrefix = 'FH2_PTA_12_8_21_';
% OutputMaskPrefix = 'Mask_';
% OutputTissuePrefix = 'Tissue_';
% Njm = 2048; Nim = 2048; Nkm = 1024;
% MRange = 141+[100:1123];
% NewSize = [512,512,256];

% Read in mask data and tissue data
fprintf('... reading images ...\n');
kstart = MRange(1); kend = MRange(Nkm); 
I3DM = uint8(false(Njm,Nim,Nkm)); 
% I3DT = uint8(false(Njm,Nim,Nkm)); 
for k=kstart:kend
  fnamein = sprintf('%s%s%04d.png',InputMaskPath,InputMaskPrefix,k);
  I3DM(:,:,k-kstart+1) = 255*uint8(imread(fnamein));
%   fnamein = sprintf('%s%s%04d.png',InputTissuePath,InputTissuePrefix,k);
%   I3DT(:,:,k-kstart+1) = uint8(imread(fnamein));
end

% resize
fprintf('... eroding image ...\n');
I3DM = imerode(I3DM,ones(3,3,3));
% I3DMCut =  I3DM;
I3DMCut2 = I3DM;
for k=1:size(I3DMCut2,3)
%     I3DMCut(:,:,k) = tril(I3DMCut(:,:,k));
    I3DMCut2(:,:,k) = triu(I3DMCut2(:,:,k));
end

fprintf('... resizing images ...\n');
% IMR = imresize3(I3DM,NewSize,'Antialiasing',true);
% IMCR = imresize3(I3DMCut,NewSize,'Antialiasing',true);
IMCR2 = imresize3(I3DMCut2,NewSize,'Antialiasing',true);
% ITR = imresize3(I3DT,NewSize,'Antialiasing',true);

% write
fprintf('... writing images ...\n');
parfor k=1:NewSize(3)
%      fnameout = sprintf('%s%s%04d.png',OutputMaskPath,OutputMaskPrefix,k);
%      imwrite(IMR(:,:,k),fnameout);
%      fnameout = sprintf('%s%s%04d.png',OutputMaskPath2,OutputMaskPrefix,k);
%      imwrite(IMCR(:,:,k),fnameout);
     fnameout = sprintf('%s%s%04d.png',OutputMaskPath3,OutputMaskPrefix,k);
     imwrite(IMCR2(:,:,k),fnameout);
%      fnameout = sprintf('%s%s%04d.png',OutputTissuePath,OutputTissuePrefix,k);
%      imwrite(ITR(:,:,k),fnameout);
end
