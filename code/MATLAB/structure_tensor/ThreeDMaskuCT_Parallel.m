%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ThreeDMaskuCT_Parallel.m
%
% This script loads in uCT images that have been resliced orthogonal to 
% the LV chamber long axis. This means the long axis of the LV chamber is 
% vertical and the LV and RV chamber centroids lie along a horizontal line
% when the image slices are displayed.
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
% Lower memory requirements and faster processing over all is obtained by 
% reading and processing strips of image slices. To reduce memory overhead
% there are few global variables.
%
% To clean up the masks (i.e. remove unconnected parts in 3D) the full 
% 3D image is required. To save memory both 
%
% Modified by: Mark Trew, November 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

t0 = cputime;
tStart = tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% User defined parameters.
%
% It would be nice to have one place for picking up the user defined
% parameters for each problem. On the other hand, having them appear as
% comments in the actual script file does help with documenting where data
% was generated from and with what settings.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Anandita example - update this.
% H1C2H dry
kindex = [1:1420];
jindex = [1:1596];
iindex = [1:1956]; % truncate on the right of the image
FileTemplate = 'H1C2H_';
InputFileExtension = 'png';
ImageFolder = '../test4/Images/';
TopLevelMaskFolder = '../test4/Masks/'; % set to empty if no output desired
TopLevelDisplayFolder = '../test4/Display/'; % set to empty if no output desired
VentricleMaskFolder = 'Ventricles/';
LVMaskFolder = 'LVCavity/';
RVMaskFolder = 'RVCavity/';
MaskedImageFolder = 'MaskedTissue/';

Nj = length(jindex); Ni = length(iindex); Nk = length(kindex);

% Imaging type and processing parameters 
ImageClass = 'uint8'; MaxVal = 2^8;
DigitsInImageSequence = 4; % number of digits in image numbering pattern
TopPadding = 20; % pad top of image stack out by 20 blank layers.
BottomPadding = 20; % bottom of image stack also padded by 20 blank layers.
MedianFilterWidth = 5;
ClosingWidth = 5; % use 5 or airdry or 11 for PTA images
OpeningWidth = 11;
GFiltSD = 7; % Gaussian filter standard deviation
Threshold = 100;

% Set up for parallel processing. For 1440 slices, 
Overlap = 25; % after some experimentation, this seems adequate.
% choose stride size so that no more that 12 processes are used. Use 
% ceiling so that stride is at least this number.
NThreads = 12;
Stride = ceil(Nk/NThreads);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parallel loop. Problem is broken up into chunks and each is processed independently.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for i=1:ceil(Nk/Stride)
parfor i=1:ceil(Nk/Stride)
    
    % Parallel range and new NkP 
    ks = max((i-1)*Stride-Overlap+1,1);
    ke = min(i*Stride+Overlap,Nk);
    NkP = length(ks:ke);
    kindexP = [ks:ke];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Set up and load in image stack
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Allocate memory
    fprintf('... allocating memory strip %d ...\n',i);
    I = cast(zeros(Nj,Ni,NkP),ImageClass);

    % Load in image set
    fprintf('... loading images strip %d ...\n',i);
    if not(isfolder(ImageFolder)) fprintf('Folder: %f does not exist.\n',ImageFolder); end
    fstring = sprintf('%s%s%%0%dd.png',ImageFolder,FileTemplate,DigitsInImageSequence);
    for k=1:NkP
        fnamein = sprintf(fstring,kindexP(k));
        S = imread(fnamein);
        I(:,:,k) = S(jindex,iindex);
    end

    % Add on zero padding to top of image beyond truncation if it is the
    % last strip of images
    if i == ceil(Nk/Stride)
        I = padarray(I,[0,0,TopPadding],0,'post');
        kindexP = [kindexP,max(kindexP)+[1:TopPadding]];
        NkP = length(kindexP);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Extract tissue mask
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Adjust intensity signals, elevate from background, morphological
    % operate, filter and threshold to find a mask. This simple approach
    % works due to consistency of uCT imaging through choice of imaging
    % parameters to match a common intensity histogram.
    fprintf(sprintf('... adjusting and saturating intensity - strip %d...\n',i));
    Mask = uint8(single(imadjustn(I)).^1.2);
    fprintf(sprintf('... median filter  - strip %d...\n',i));
    Mask = medfilt3(Mask,[MedianFilterWidth,MedianFilterWidth,MedianFilterWidth]);
    fprintf(sprintf('... morphological operations  - strip %d...\n',i));
    Mask = imopen(imclose(Mask,ones(ClosingWidth,ClosingWidth,ClosingWidth)),ones(OpeningWidth,OpeningWidth,OpeningWidth));
    fprintf(sprintf('... Gaussian filtering  - strip %d...\n',i));
    Mask = imgaussfilt3(Mask,GFiltSD);
    fprintf(sprintf('... Thresholding  - strip %d...\n',i));
    Mask = Mask > Threshold;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Identify cavities
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(sprintf('... Masking cavity  - strip %d...\n',i));
    CavityMask = Mask;
    for k=1:NkP
        CavityMask(:,:,k) = logical(imfill(CavityMask(:,:,k),'holes'));
    end
    CavityMask = logical(imsubtract(CavityMask,Mask));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Define output range for strip
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i==1
        StartLayer = 1;
        EndLayer = min(Stride,NkP);
    elseif i==ceil(Nk/Stride)
        StartLayer = Overlap+1;
        EndLayer = NkP;
    else
        StartLayer = Overlap+1;
        EndLayer = Stride+Overlap;
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Write mask strips
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(TopLevelMaskFolder)
        fprintf('... writing mask, cavity and masked tissue image strip %d ...\n',i);
        if not(isfolder(TopLevelMaskFolder)) mkdir(TopLevelMaskFolder); end
       
        % Write mask
        fprintf('    1. writing mask image strip %d ...\n',i);    
        if not(isfolder([TopLevelMaskFolder,VentricleMaskFolder])) mkdir([TopLevelMaskFolder,VentricleMaskFolder]); end
        fstring = sprintf('%s%s%%0%dd.png',[TopLevelMaskFolder,VentricleMaskFolder],FileTemplate,DigitsInImageSequence);
        for k=StartLayer:EndLayer
            fnameout = sprintf(fstring,kindexP(k)); 
            imwrite(Mask(:,:,k),fnameout,'BitDepth',1);
        end

        % Write cavity mask - this is temporarily written to the LV folder
        fprintf('    2. writing cavity mask image strip %d ...\n',i);
        if not(isfolder([TopLevelMaskFolder,LVMaskFolder])) mkdir([TopLevelMaskFolder,LVMaskFolder]); end
        fstring = sprintf('%s%s%%0%dd.png',[TopLevelMaskFolder,LVMaskFolder],FileTemplate,DigitsInImageSequence);
        for k=StartLayer:EndLayer
            fnameout = sprintf(fstring,kindexP(k)); 
            imwrite(CavityMask(:,:,k),fnameout,'BitDepth',1);
        end
    end
       
end % end the parallel loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read masks back in and clean up. This is retaining the single largest 
% mask in the ventricles and the two largest masks in the cavity. Also
% write the two largest masks separately.
%
% It might be best not to write the masked tissue above and work on
% cleaning the mask first.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(TopLevelMaskFolder)
    fprintf('... allocating image and mask memory ...\n',i);
    I = cast(zeros(Nj,Ni,Nk+TopPadding),ImageClass);
    M = cast(zeros(Nj,Ni,Nk+TopPadding),'logical');
    LVM = cast(zeros(Nj,Ni,Nk+TopPadding),'logical');
    RVM = cast(zeros(Nj,Ni,Nk+TopPadding),'logical');

    % Load in original image set and trim as specified in jindex and iindex
    % Load in masks
    fprintf('... loading images and masks ...\n',i);
    fstringImage = sprintf('%s%s%%0%dd.png',ImageFolder,FileTemplate,DigitsInImageSequence);
    fstringMask = sprintf('%s%s%%0%dd.png',[TopLevelMaskFolder,VentricleMaskFolder],FileTemplate,DigitsInImageSequence);
    fstringLVMask = sprintf('%s%s%%0%dd.png',[TopLevelMaskFolder,LVMaskFolder],FileTemplate,DigitsInImageSequence);
    parfor k=1:Nk
        fnamein = sprintf(fstringImage,kindex(k));
        S = imread(fnamein);
        I(:,:,k) = S(jindex,iindex);
        fnamein = sprintf(fstringMask,kindex(k));
        M(:,:,k) = imread(fnamein);
        fnamein = sprintf(fstringLVMask,kindex(k));
        LVM(:,:,k) = imread(fnamein);
    end

    % Clean ventricle mask
    fprintf('... cleaning ventricle mask ...\n',i);
    stats = regionprops(M,'Area','PixelIdxList');
    [Order,Oloc]=sort(cat(1,stats(:).Area),'descend');
    M(:) = false;
    M(stats(Oloc(1)).PixelIdxList) = 1;
    
    % Clean cavity mask
    fprintf('... cleaning LV and RV cavity mask ...\n',i);
    stats = regionprops(LVM,'Area','PixelIdxList');
    [Order,Oloc]=sort(cat(1,stats(:).Area),'descend');
    LVM(:) = false;
    LVM(stats(Oloc(1)).PixelIdxList) = 1;
    RVM(stats(Oloc(2)).PixelIdxList) = 1;

    % Rewrite the files 
    fprintf('... Write cleaned masks and masked images...\n',i);

    if not(isfolder([TopLevelMaskFolder,VentricleMaskFolder])) mkdir([TopLevelMaskFolder,VentricleMaskFolder]); end
    if not(isfolder([TopLevelMaskFolder,LVMaskFolder])) mkdir([TopLevelMaskFolder,LVMaskFolder]); end
    if not(isfolder([TopLevelMaskFolder,RVMaskFolder])) mkdir([TopLevelMaskFolder,RVMaskFolder]); end
    if not(isfolder([TopLevelMaskFolder,MaskedImageFolder])) mkdir([TopLevelMaskFolder,MaskedImageFolder]); end

    fstringMask = sprintf('%s%s%%0%dd.png',[TopLevelMaskFolder,VentricleMaskFolder],FileTemplate,DigitsInImageSequence);
    fstringLVMask = sprintf('%s%s%%0%dd.png',[TopLevelMaskFolder,LVMaskFolder],FileTemplate,DigitsInImageSequence);
    fstringRVMask = sprintf('%s%s%%0%dd.png',[TopLevelMaskFolder,RVMaskFolder],FileTemplate,DigitsInImageSequence);
    fstringMaskedTissue = sprintf('%s%s%%0%dd.png',[TopLevelMaskFolder,MaskedImageFolder],FileTemplate,DigitsInImageSequence);
    parfor k=1:Nk+TopPadding
        fnameout = sprintf(fstringMask,k); 
        imwrite(M(:,:,k),fnameout);
        fnameout = sprintf(fstringLVMask,k); 
        imwrite(LVM(:,:,k),fnameout);
        fnameout = sprintf(fstringRVMask,k); 
        imwrite(RVM(:,:,k),fnameout);
        fnameout = sprintf(fstringMaskedTissue,k); 
        imwrite(uint8(M(:,:,k)).*I(:,:,k),fnameout);
    end
           
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Quantify mask features.
%
% Mask quantification includes:
% 1. Surface areas 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(TopLevelDisplayFolder)
    if not(isfolder(TopLevelDisplayFolder)) mkdir(TopLevelDisplayFolder); end
    fprintf('... quantify masks ...\n',i);

    % Total LV volume and surface area - just use 1:Nk so that the top
    % surface is not eroded and is not included in the surface area
    % calculation
    LVVolume = sum(LVM,'all');
    SALV = LVM(:,:,1:Nk) - imerode(LVM(:,:,1:Nk),ones(3,3,3));
    LVSurfaceArea = sum(SALV,'all');

    % Total RV volume and surface area
    RVVolume = sum(RVM,'all');
    SARV = RVM(:,:,1:Nk) - imerode(RVM(:,:,1:Nk),ones(3,3,3));
    RVSurfaceArea = sum(SARV,'all');
   
    % Total Ventricular volume and surface area
    VTVolume = sum(M,'all');
    SAVT = M(:,:,1:Nk) - imerode(M(:,:,1:Nk),ones(3,3,3));
    VTSurfaceArea = sum(SAVT,'all');
    
    % Slice volume and surface area
    WorkingV = reshape(LVM,[Nj*Ni,Nk+BottomPadding]); WorkingV = sum(WorkingV,1);
    WorkingA = reshape(SALV,[Nj*Ni,Nk]); WorkingA = sum(WorkingA,1);
    SliceLVStats = [[1:Nk-BottomPadding+1]',WorkingV(BottomPadding:Nk)',WorkingA(BottomPadding:Nk)'];

    WorkingV = reshape(RVM,[Nj*Ni,Nk+BottomPadding]); WorkingV = sum(WorkingV,1);
    WorkingA = reshape(SARV,[Nj*Ni,Nk]); WorkingA = sum(WorkingA,1);
    SliceRVStats = [[1:Nk-BottomPadding+1]',WorkingV(BottomPadding:Nk)',WorkingA(BottomPadding:Nk)'];

    WorkingV = reshape(M,[Nj*Ni,Nk+BottomPadding]); WorkingV = sum(WorkingV,1);
    WorkingA = reshape(SAVT,[Nj*Ni,Nk]); WorkingA = sum(WorkingA,1);
    SliceVTStats = [[1:Nk-BottomPadding+1]',WorkingV(BottomPadding:Nk)',WorkingA(BottomPadding:Nk)'];
    
    % Write data to file for later plotting
    save(sprintf('%sMaskQuantification.mat',TopLevelDisplayFolder),'LVVolume','LVSurfaceArea',...
        'RVVolume','RVSurfaceArea','VTVolume','VTSurfaceArea','SliceLVStats','SliceRVStats','SliceVTStats');
    
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Down sample the masks for display purposes. The masks are converted to
% uint8 data and eroded by 1 voxel. The new mask sizes have the first
% dimension as 512 and other dimensions are relative to this and maintain
% the original aspect ratio. Resizing is done with antialiasing so that the
% rendered iso surfaces are smooth.
%
% Exnode file for the receptacle box are also generated.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(TopLevelDisplayFolder)
    fprintf('... downsample masks and images for display ...\n',i);
    
    % Determine new image size
    JTarget = 512;
    Alpha = JTarget/Nj;
    NewSize = [JTarget,ceil(Ni*Alpha),ceil((Nk+TopPadding)*Alpha)];
    
    % Resize image masked with non-eroded ventricles
    I = imresize3(I.*uint8(M),NewSize);    
    
    % Erode masks by 1 voxel
    M = imerode(M,ones(3,3,3));
    LVM = imerode(LVM,ones(3,3,3));
    RVM = imerode(RVM,ones(3,3,3));
  
    % Resized remainder of masks converted to 8 bit
    M = imresize3(255*uint8(M),NewSize,'Antialiasing',true);
    LVM = imresize3(255*uint8(LVM),NewSize,'Antialiasing',true);
    RVM = imresize3(255*uint8(RVM),NewSize,'Antialiasing',true);
    
    % Write resized masks for display
    fprintf('... write downsample masks and images...\n',i);
    if not(isfolder(TopLevelDisplayFolder)) mkdir(TopLevelDisplayFolder); end
    
    fprintf('    1. Ventricular mask ...\n',i);    
    if not(isfolder([TopLevelDisplayFolder,VentricleMaskFolder])) mkdir([TopLevelDisplayFolder,VentricleMaskFolder]); end
    fstring = sprintf('%s%s%%0%dd.png',[TopLevelDisplayFolder,VentricleMaskFolder],FileTemplate,DigitsInImageSequence);
    for k=1:NewSize(3)
        fnameout = sprintf(fstring,k); 
        imwrite(M(:,:,k),fnameout);
    end

    fprintf('    2. LV cavity mask ...\n',i);
    if not(isfolder([TopLevelDisplayFolder,LVMaskFolder])) mkdir([TopLevelDisplayFolder,LVMaskFolder]); end
    fstring = sprintf('%s%s%%0%dd.png',[TopLevelDisplayFolder,LVMaskFolder],FileTemplate,DigitsInImageSequence);
    for k=1:NewSize(3)
        fnameout = sprintf(fstring,k); 
        imwrite(LVM(:,:,k),fnameout);
    end

    fprintf('    3. RV cavity mask ...\n',i);
    if not(isfolder([TopLevelDisplayFolder,RVMaskFolder])) mkdir([TopLevelDisplayFolder,RVMaskFolder]); end
    fstring = sprintf('%s%s%%0%dd.png',[TopLevelDisplayFolder,RVMaskFolder],FileTemplate,DigitsInImageSequence);
    for k=1:NewSize(3)
        fnameout = sprintf(fstring,k); 
        imwrite(RVM(:,:,k),fnameout);
    end

    fprintf('    4. Masked tissue ...\n',i);
    if not(isfolder([TopLevelDisplayFolder,MaskedImageFolder])) mkdir([TopLevelDisplayFolder,MaskedImageFolder]); end
    fstring = sprintf('%s%s%%0%dd.png',[TopLevelDisplayFolder,MaskedImageFolder],FileTemplate,DigitsInImageSequence);
    for k=1:NewSize(3)
        fnameout = sprintf(fstring,k); 
        imwrite(I(:,:,k),fnameout);
    end

    % Set up display box in exnode and exelem format - write to
    % TopLevelDisplayFolder. Note the ij order are reversed for the cmgui
    % image coordinates.
    fprintf('... write exnode exelem display files...\n',i);
    Vertices = [1,1,1;NewSize(2),1,1;1,NewSize(1),1;NewSize(2),NewSize(1),1;1,1,NewSize(3);NewSize(2),1,NewSize(3);1,NewSize(1),NewSize(3);NewSize(2),NewSize(1),NewSize(3)];
    WriteSingleVolumeElement(Vertices,TopLevelDisplayFolder);
end

t1 = cputime;
tEnd = toc(tStart);
fprintf(' Total time for processing: %f s (%f min) \n',tEnd,tEnd/60);
