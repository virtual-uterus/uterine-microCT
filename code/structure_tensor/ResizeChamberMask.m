% This script loads in chambers mask, resizes and writes as 8 bit resized
% for display

% Chamber mask images
% InputPath = '../Images/ChamberMask_2048x2048x1024/';
% OutputPath = '../Images/ChamberMask_8bit_512x512x256/';
% InputPrefix = 'ChamberMask_';
% OutputPrefix = 'ChamberMask_';
% Njm = 2048; Nim = 2048; Nkm = 1024;
% MRange = 141+[100:1123];

% Chamber mask images - upper ventricles
InputPath = '../Images/Mask_Part2_CleanMask/ChamberMask_2048x2048x512/';
OutputPath = '../Images/Mask_Part2_CleanMask/ChamberMask_8bit_512x512x128/';
InputPrefix = 'ChamberMask_';
OutputPrefix = 'ChamberMask_';
Njm = 2048; Nim = 2048; Nkm = 512;
MRange = [1150:1661];

% Read in mask data, resize and writ
fprintf('... reading images ...\n');
kstart = MRange(1); kend = MRange(Nkm); 
I3D = uint8(false(Njm,Nim,Nkm)); 
for k=kstart:kend
  fnamein = sprintf('%s%s%04d.png',InputPath,InputPrefix,k);
  I3D(:,:,k-kstart+1) = 255*uint8(imread(fnamein));
end

% resize
fprintf('... resizing images ...\n');
%IR = imresize3(I3D,[512,512,256],'Antialiasing',true);
IR = imresize3(I3D,[512,512,128],'Antialiasing',true);

% write
fprintf('... writing images ...\n');
%parfor k=1:256
parfor k=1:128
     fnameout = sprintf('%s%s%04d.png',OutputPath,OutputPrefix,k);
     imwrite(IR(:,:,k),fnameout);
end
