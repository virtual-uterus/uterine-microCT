%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DiffusionTissueExtrapolation.m
%
% This script loads in a masked tissue image stack and identifies the
% tissue as having values > threshold. The tissue mask is heavily
% eroded to ensure that tissue intensity values provide boundary
% conditions for diffusion extrapolation into non-tissue regions.
%
% The diffusion extrapolation iterates until the mean of the differences
% in intensity between voxels a specified distance from the tissue
% and the nearest tissue values crosses zero.
%
% Written by: Mark Trew, November 2011.
% Modified by: Mark Trew, March 2021
% Modified by: Mark Trew, September 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% User defined parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image class, dimension and numbering sequence
%ImageClass = 'logical';
%ImageClass = 'uint16'; MaxVal = 2^16;
ImageClass = 'uint8'; MaxVal = 2^8;
% Ni = 2048; Nj = 2048; Nk = 2048; 
% kindex = [141:2188]; % 
% 
% % Image input and output paths
% ImageInput = '../Images/Masked/';
% InputFileTemplate = 'FH2_PTA_12_8_21_';
% ImageOutput = '../Images/MaskedExtrapolated/';
% OutputFileTemplate = 'FH2_PTA_12_8_21_';
% DigitsInImageSequence = 4;

% Upper part of heart - slices 1150 to 1661, i.e. 512 images
% Ni = 2048; Nj = 2048; Nk = 512; 
% kindex = [1150:1661]; % 
% ImageInput = '../Images/Mask_Part2_CleanMask/AllFramesMasked/';
% InputFileTemplate = 'FH2_PTA_12_8_21_';
% ImageOutput = '../Images/Mask_Part2_CleanMask/AllFramesMasked_Extrapolated/';
% OutputFileTemplate = 'FH2_PTA_12_8_21_';
% DigitsInImageSequence = 4;

% Heart 1
% Nj = 2048; Ni = 2352; 
% kindex = [260:1752];
% Nk = length(kindex);
% InputFileTemplate = 'FH1_PTA_20_7_21_';
% ImageInput = '../ImagesHeart1/MaskedEdited/';
% ImageOutput = '../ImagesHeart1/MaskedEdited_Extrapolated/';
% OutputFileTemplate = 'FH1_PTA_20_7_21_';
% DigitsInImageSequence = 4;

% Heart 2 Air Dry
% Nj = 1237; Ni = 1237;
% kindex = [548:1060];
% Nk = length(kindex);
% InputFileTemplate = 'AirDry_';
% ImageInput = '../ImagesHeartDry/Masked/';
% ImageOutput = '../ImagesHeartDry/MaskedExtrapolated/';
% OutputFileTemplate = 'AirDry_';
% DigitsInImageSequence = 4;

% Images from Air dry H1C1H
Nj = 1589; Ni = 1771;
kindex = [370:1689];
Nk = length(kindex);
ImageInput = '../Images_H1C1H/MaskedClean/';
InputFileTemplate = 'H1C1H_';
ImageOutput = '../Images_H1C1H/MaskedCleanExtrapolated/';
OutputFileTemplate = 'H1C1H_';
DigitsInImageSequence = 4;

% Tissue mask erosion threshold and radius
TissueMaskErosionThreshold = 20;
TissueMaskErosionRadius = 4;
%PadWidth = 2;
PadWidth = 0;

% Set the tissue boundary diffusion testing distance (voxels) - points 
% at this distance from the tissue will be used to assess the stopping 
% criteria for the diffusion iterations.
TissueBoundaryDiffusionTestingDistance = 7; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set up and load in image stack
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory
fprintf('... allocating memory ...\n');
I = cast(zeros(Nj,Ni,Nk),ImageClass);

% Load in image set
fprintf('... loading images ...\n');
fstring = sprintf('%s%s%%0%dd.png',ImageInput,InputFileTemplate,...
    DigitsInImageSequence);
parfor k=1:length(kindex)
  if ~mod(k,100), fprintf(' image: %d\n',k); end
  fnamein = sprintf(fstring,kindex(k)); 
  I(:,:,k) = imread(fnamein);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Identify tissue mask and erode.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('...mask and erosion...\n');
% 1. Identify tissue mask
IB = (I >= TissueMaskErosionThreshold);
% 2. erode tissue mask - the erosion filter is spherical with diameter
%    2*R+1 voxels
R = TissueMaskErosionRadius;
NHood = zeros(2*R+1,2*R+1,2*R+1);
NHood(R+1,R+1,R+1) = 1;
D = bwdist(NHood);
NHood = (D <= R);
IB = imerode(IB,NHood);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Find voxels a specified distance from the tissue. Store nearest tissue 
% values.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('...finding Euclidian distances...\n');
DT = TissueBoundaryDiffusionTestingDistance;
[D,L]=bwdist(IB);
TIdx = find(D >= DT-0.5 & D <= DT+0.5);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set up and perform diffusion.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set timing initialization
trun0 = cputime;

% Identify the points that will be used as boundary conditions 
% on the diffusion 
FixedPoints = find(IB);

% "time" step
Delta = 0.3;

% Set up a 3D 7 point diffusion filter and then construct a circulant 
% matrix. To fit with previous use of this approach, a 5x5x5 template is used.
F = zeros(5,5,5);
F(2,2,2) = -6*Delta + 1; % include contribution from identity matrix
F(2-1,2,2) = 1*Delta; F(2+1,2,2) = 1*Delta;
F(2,2-1,2) = 1*Delta; F(2,2+1,2) = 1*Delta;
F(2,2,2-1) = 1*Delta; F(2,2,2+1) = 1*Delta;

% Circular shift offset - this is the distance to the middle
% of the 5x5x5 template
CS = -(2*(Ni+2*PadWidth)*(Nj+2*PadWidth)+2*(Ni+2*PadWidth)+3);

% Full diffusion weight matrix
WD= single(zeros([Ni+2*PadWidth,Nj+2*PadWidth,Nk+2*PadWidth]));
WD(1:5,1:5,1:5) = F;
WD = reshape(WD,[(Ni+2*PadWidth)*(Nj+2*PadWidth)*(Nk+2*PadWidth),1]);
WD = flipud(circshift(WD,CS));
WDf = single(fft(WD)); clear WD;

% Initial conditions on diffused image DI. Delta is the effective 
% "time step". If it is too large the diffusion will be unstable.
% Clear the I and IB arrays as they are no longer needed.
DI = I.*cast(IB,'uint8'); clear I IB D;
DI = single(reshape(DI,Ni*Nj*Nk,1));
BVal = DI(FixedPoints);

% Iterate over diffusion steps while the sign of the mean difference 
% between the diffused image intensity DT voxels from the tissue and the 
% nearest boundary condition is less than zero. The L array is the index of 
% the original closest non-zero pixel to the test group, i.e. the nearest 
% boundary condition.
% Note: this algorithm could be rewritten to potentially be faster. 
% Using the same tools as for structure tensor analysis the F could be 
% set up as a circulant matrix. Then the update operation would be:
% DI = Delta*F*DI + DI = (Delta*F+I)*DI = S*DI
% Here S becomes a new circulant matrix and the S*DI operation simply 
% becomes a convolution operation and should be relatively fast.
n=1;
md = -1;
while  md <= 0
  md = mean(double(DI(TIdx))-double(DI(L(TIdx))));
  if ~mod(n,1), fprintf('... n: %d, mean diff: %f\n',n,md); end
  
  % apply diffusion filter to current state of diffused image. 
  DIf = single(fft(DI));   
  DI = ifft(DIf .* WDf);
  % apply fixed boundary conditions
  DI(FixedPoints) = BVal;
  % update counter
  n = n+1;
end

trun1 = cputime;
fprintf(' Total time for %d its: %f\n',n,trun1-trun0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Write boundary smoothed image files.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('... writing images ...\n');
fstring = sprintf('%s%s%%0%dd.png',ImageOutput,OutputFileTemplate,DigitsInImageSequence);
for k=1:length(kindex),
  if ~mod(k,100), fprintf(' image: %d\n',k); end;
  fnameout = sprintf(fstring,kindex(k)); 
  imwrite(DI(:,:,k),fnameout,'png');
end;
