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
DigitsInImageSequence = 3;

%%
kindex = 80:534; 
Nj = 548; Ni = 1124;
Nk = length(kindex);
base_dir = join([getenv("HOME"), "Documents/phd"], '/');
ImageInput = base_dir + '/AWA015_PTA_1_Rec_Trans/downsampled/ST/masked/';
InputFileTemplate = 'AWA015_PTA_1_';
ImageOutput = base_dir + '/AWA015_PTA_1_Rec_Trans/downsampled/ST/extrapolated/';
OutputFileTemplate = 'AWA015_PTA_1_';
img_digit_size = '%s%s%03d.%s';
%%
% kindex = 1:325; 
% Nj = 518; Ni = 753;
% Nk = length(kindex);
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% ImageInput = base_dir + '/AWA015_PTA_2_Ova_Rec_Trans/downsampled/ST/masked/';
% InputFileTemplate = 'AWA015_PTA_2_Ova_Trans_';
% ImageOutput = base_dir + '/AWA015_PTA_2_Ova_Rec_Trans/downsampled/ST/extrapolated/';
% OutputFileTemplate = 'AWA015_PTA_2_Ova_Trans_';
% img_digit_size = '%s%s%03d.%s';
%% Downsampled version
% Nj = 469; Ni = 512;
% kindex = 1:138;
% Nk = length(kindex);
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% ImageInput = base_dir + '/AWA014_PTA_2_Rec_Trans/downsampled/ST/masked/';
% InputFileTemplate = 'AWA014_PTA_2_';
% ImageOutput = base_dir + '/AWA014_PTA_2_Rec_Trans/downsampled/ST/extrapolated/';
% OutputFileTemplate = 'AWA014_PTA_2_';

%% Full resolution segment version
% kindex = 830:1030; 
% Nj = 1401; Ni = 1001;
% Nk = length(kindex);
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% ImageInput = base_dir + '/AWA014_PTA_2_Rec_Trans/ST/masked/';
% InputFileTemplate = 'AWA014_PTA_2_';
% ImageOutput = base_dir + '/AWA014_PTA_2_Rec_Trans/ST/extrapolated/';
% OutputFileTemplate = 'AWA014_PTA_2_';


%% Reoriented downsampled version
% Nj = 367; Ni = 473;
% kindex = 1:352;
% Nk = length(kindex);
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% ImageInput = base_dir + '/AWA014_PTA_2_Rec_Trans/downsampled/ST/masked/';
% InputFileTemplate = 'AWA014_PTA_2_';
% ImageOutput = base_dir + '/AWA014_PTA_2_Rec_Trans/downsampled/ST/extrapolated/';
% OutputFileTemplate = 'AWA014_PTA_2_';


%% Reoriented full resolution segment version
% Nj = 420; Ni = 388;
% kindex = 1:383;
% Nk = length(kindex);
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% ImageInput = base_dir + '/AWA014_PTA_2_Rec_Trans/ST/masked/';
% InputFileTemplate = 'AWA014_PTA_2_';
% ImageOutput = base_dir + '/AWA014_PTA_2_Rec_Trans/ST/extrapolated/';
% OutputFileTemplate = 'AWA014_PTA_2_';


% Tissue mask erosion threshold and radius
TissueMaskErosionThreshold = 10;
TissueMaskErosionRadius = 2;

% Set the tissue boundary diffusion testing distance (voxels) - points 
% at this distance from the tissue will be used to assess the stopping 
% criteria for the diffusion iterations.
TissueBoundaryDiffusionTestingDistance = 2; 

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

% 3D 7 point diffusion filter
F = zeros(3,3,3);
F(2,2,1)=1; 
F(2,1,2)=1; F(1,2,2)=1; F(3,2,2)=1; F(2,3,2)=1; F(2,2,2)=-6; 
F(2,2,3)=1;

% Set up a 3D 7 point diffusion filter and then construct a circulant 
% matrix.

% Initial conditions on diffused image DI. Delta is the effective 
% "time step". If it is too large the diffusion will be unstable.
% Clear the I and IB arrays as they are no longer needed.
DI = I.*cast(IB,'uint8'); clear I  IB  D;
Delta = 0.3;

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
  if ~mod(n,10), fprintf('... n: %d, mean diff: %f\n',n,md); end
  % apply diffusion filter to current state of diffused image. DI is 
  % exterior padded with replicated values to enforce zero flux b.c.
  DIhat = imfilter(padarray(DI,[1,1,1],'replicate'),F);
  DIhat = DIhat(2:Nj+1,2:Ni+1,2:Nk+1);
  % apply fixed boundary conditions
  DIhat(FixedPoints) = 0;
  % update the diffused image
  DI = Delta*DIhat + DI;
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
for k=1:length(kindex)
  if ~mod(k,100), fprintf(' image: %d\n',k); end
  fnameout = sprintf(fstring,kindex(k)); 
  imwrite(DI(:,:,k),fnameout,'png');
end
