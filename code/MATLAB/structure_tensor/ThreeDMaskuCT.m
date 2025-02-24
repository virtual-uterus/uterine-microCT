%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ThreeDMaskuCT.m
%
% This script loads in uCT images and then masks them.
%
% It is expected that the uCT images will have first been placed in a
% consistent orientation where the long axis of the LV is vertical and the
% LV and RV chamber centroids lie along a horizontal line when the image
% slices are displayed.
%
% The relatively simple approach to segmentation is possible because the
% uCT imaging settings are altered for each sample to match a common
% histogram profile.
%
% For consistency between hearts, not only are the ventricles reoriented
% along the LV axis, but the image stack is truncated so that a height of
% 1.4 mm from the apex of the heart is retained. Valve planes and atria are
% not retained for later comparison and analysis.
%
% Updated by: 
%    Mark Trew, November 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

% Load 8 bit uCT images
% process and filter and threshold to generate mask.
ImageClass = 'uint8'; MaxVal = 2^8;
DigitsInImageSequence = 3; % number of digits in image numbering pattern
TopPadding = 0; % pad top of image stack out by 20 blank layers.

%% Downsampled version
% kindex = 1:138;
% iindex = 1:512;
% jindex = 1:469; % truncate on the right of the image
% ImageInput = '../AWA014_PTA_2_Rec_Trans/downsampled/';
% InputFileTemplate = 'AWA014_PTA_2_';
% ImageOutput = ImageInput;
% OutputTemplate = 'AWA014_PTA_2_';
% InputFileExtension = 'png';
% Nj = length(jindex); Ni = length(iindex); Nk = length(kindex);

%% Code edits
% Changed the overlap and stride to match the one from the reslice code.
%%
kindex = 1:325; 
iindex = 1:753;
jindex = 1:518;
base_dir = join([getenv("HOME"), "Documents/phd"], '/');
ImageInput = base_dir + '/AWA015_PTA_2_Ova_Rec_Trans/downsampled/';
InputFileTemplate = 'AWA015_PTA_2_Ova_Trans_';
ImageOutput = ImageInput + 'ST/';
OutputTemplate = 'AWA015_PTA_2_Ova_Trans_';
InputFileExtension = 'png';
Nj = length(jindex); Ni = length(iindex); Nk = length(kindex);
img_digit_size = '%s%s%03d.%s';
%%
% kindex = 1:256; 
% iindex = 1:512;
% jindex = 1:250;
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% ImageInput = base_dir + '/AWA015_PTA_1_Rec_Trans/downsampled/';
% InputFileTemplate = 'AWA015_PTA_1_';
% ImageOutput = ImageInput + 'ST/';
% OutputTemplate = 'AWA015_PTA_1_';
% InputFileExtension = 'png';
% Nj = length(jindex); Ni = length(iindex); Nk = length(kindex);
% img_digit_size = '%s%s%03d.%s';
%% Full resolution segment version
% kindex = 830:1030; 
% iindex = 220:1220;
% jindex = 50:1450;
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% ImageInput = base_dir + '/AWA014_PTA_2_Rec_Trans/';
% InputFileTemplate = 'AWA014_PTA_2_';
% ImageOutput = ImageInput + 'ST/';
% OutputTemplate = 'AWA014_PTA_2_';
% InputFileExtension = 'bmp';
% Nj = length(jindex); Ni = length(iindex); Nk = length(kindex);
% img_digit_size = '%s%s%08d.%s';

%% Reoriented full resolution segment version
% kindex = 1:1243;
% jindex = 1:482;
% iindex = 1:1043; % truncate on the right of the image
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% ImageInput = base_dir + '/AWA014_PTA_2_Rec_Trans/ST/';
% InputFileTemplate = 'AWA014_PTA_2_';
% ImageOutput = ImageInput;
% OutputTemplate = 'AWA014_PTA_2_';
% InputFileExtension = 'png';
% Nj = length(jindex); Ni = length(iindex); Nk = length(kindex);
% img_digit_size = '%s%s%%0%dd.%s';

%% Reoriented downsampled version
% kindex = 1:256;
% jindex = 1:;
% iindex = 1:56; % truncate on the right of the image
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% ImageInput = base_dir + '/AWA014_PTA_2_Rec_Trans/downsampled/';
% InputFileTemplate = 'AWA014_PTA_2_';
% ImageOutput = ImageInput + "ST/";
% OutputTemplate = 'AWA014_PTA_2_';
% InputFileExtension = 'png';
% Nj = length(jindex); Ni = length(iindex); Nk = length(kindex);
% img_digit_size = '%s%s%%0%dd.%s';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Allocate memory and load in image stack
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory
fprintf('... allocating memory ...\n');
I = cast(false(Nj,Ni,Nk),ImageClass);

% Load in image set
fprintf('... loading images ...\n');
parfor k=1:length(kindex)
    fnamein = sprintf(img_digit_size,ImageInput,InputFileTemplate,kindex(k),InputFileExtension);
    S = imread(fnamein);
    I(:,:,k) = S(jindex,iindex);
end

% Add on zero padding to top of image beyond truncation
I = padarray(I,[0,0,TopPadding],0,'post');
kindex = kindex:(max(kindex)+TopPadding);
Nk = length(kindex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pre-process image stack and threshold for segmentation. Break image stack
% into blocks or strips and process each in parallel and then reassemble at
% the end.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for strips and working memory allocation
Overlap = 10;
NThreads = 12;
Stride = ceil(Nk/NThreads);
Working = cell(ceil(Nk/Stride),1);
% Stride = 20;
% Working = cell(ceil(Nk/Stride),1);

% Parallel process of independent strips
fprintf(sprintf('... Processing and masking ...\n',i));
parfor i=1:ceil(Nk/Stride)
    Is = max((i-1)*Stride-Overlap+1,1);
    Ie = min(i*Stride+Overlap,Nk);
    fprintf(sprintf('... adjusting and saturating intensity - split %d...\n',i));
    Working{i} = uint8(single(imadjustn(I(:,:,Is:Ie))).^1.2);
    fprintf(sprintf('... median filter  - split %d...\n',i));
    Working{i} = medfilt3(Working{i},[5,5,5]);
    fprintf(sprintf('... morphological operations  - split %d...\n',i));
    Working{i} = imopen(imclose(Working{i},ones(5,5,5)),ones(11,11,11));
%    Working = imopen(imclose(Working,ones(11,11,11)),ones(11,11,11)); % seems required for the PTA images
    fprintf(sprintf('... Gaussian filtering  - split %d...\n',i));
    Working{i} = imgaussfilt3(Working{i},7);
    fprintf(sprintf('... Thresholding  - split %d...\n',i));
    Working{i} = Working{i} > 100;
end
% Serial reassembly into the mask array
fprintf(sprintf('... Assembling strips ...\n',i));
IPFM = cast(false(Nj,Ni,Nk),'logical');
StoreEnd = min(Stride,Nk);
IPFM(:,:,1:StoreEnd) = Working{1}(:,:,1:Stride);
for i=2:ceil(Nk/Stride)-1
    StoreStart = max((i-1)*Stride+1,1);
    StoreEnd = min(i*Stride,Nk);
    IPFM(:,:,StoreStart:StoreEnd) = Working{i}(:,:,Overlap+1:end-Overlap);
end
StoreStart = max((ceil(Nk/Stride)-1)*Stride+1,1);
IPFM(:,:,StoreStart:Nk) = Working{ceil(Nk/Stride)}(:,:,Overlap+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Output masks and tissue.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('... writing images ...\n');
% Write mask
parfor k=1:Nk
    imwrite(IPFM(:,:,k),sprintf('%smask/%s%03d.png',ImageOutput,OutputTemplate,kindex(k)),'BitDepth',1);
end

% Write masked tissue
parfor k=1:Nk
    imwrite(uint8(IPFM(:,:,k)).*I(:,:,k),sprintf('%smasked/%s%03d.png',ImageOutput,OutputTemplate,kindex(k)));
end
