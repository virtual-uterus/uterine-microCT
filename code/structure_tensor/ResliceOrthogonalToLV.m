%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ResliceOrthogonalToLV.m
%
% This script loads in mask files and segments out the ventricles.
%
% LV lumen principal axis is found (vector a). The LV centroid and RV centroid are
% found. A plane orthogonal to the principal axis is defined with its +j direction 
% oriented along the line connecting the LV centroid to the LV principal
% axis from the RV centroid (vector c) and its remaining axis as b = c x a.
% 
% The background masked raw images are resliced using these directions.
%
% The resulting image stack has b and c vectors as i and j and k is along
% the a vector direction. This is also the helix plane on which the fiber
% orientations are projected for a scaler summary of the tissue structure.
%
% The resampled image plane is tight around the ventricles with a buffer of 
% 20 pixels in each direction. There are an additional 20 planes at the
% apex and base the complete the isotropic padding.

% The resampled planes and number of images are at an identical scale to
% the original images. For these fetal hearts the image plane pixels are
% 1um x 1um and the distance between planes is also 1um.
%
% The raw images are resampled using linear interpolation. An alternative
% would be nearest neighbour sampling. Given the resampling is at a similar 
% scale to the raw images, the linear interpolation should not be
% excessively smoothing and it will not have the rapid changes in intensity
% gradients associated with nearest neighbour.
%
% Updated by: 
%    Mark Trew, June 2022.
%    Mark Trew, July 2022.
%    Mark Trew, November 2022
%    Mark Trew, December 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Data and path locations
%
% If the MaskPath is left empty, a new mask is generated and this is used
% to isolate the LV cavity and find the centroids of both the LV and RV
% cavities. If MaskPath is no empty then an existing mask is used.
% Generating a new mask is much slower. However, it reduces the amount of
% intermediate information that is stored and arguably makes the pipeline
% cleaner.
%
% If ImageOutput or DataOutput are empty strings then the respective 
% resliced image stack and the axes figure are not written.
%
% Note: to set the i, j and k ranges Fiji/ImageJ is useful as a box can be
% drawn around the ventricles and the stack moved forward and back to
% ensure the tissue remains within the bounding box. (i,j,k) are in image
% coordinates with the origin in the upper left corner of the image (i.e.
% coordinates 1 below). While Matlab has its display image coordinates in
% image space (i.e. 1 below), the data access in arrays is made with the
% coordinates in 2 below. The data coordinate system is right handed.
%
%      (1)               (2) 
%      ^ k               ^ k
%     /                 /   
%    /                 /   
%   +----> i          +----> j (iindex)
%   |                 | 
%   |                 |
%   v j               v i (jindex)
%
% A new variable PermuteOrder is added. This is used to 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Full downsampled version
%% Code edits
% Changed the overlap from 25 to 10 for code to work, don't know how it
% affects things.
% Removed the permutation section of the code as well.
%%
%% Downsampled version
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% ImagePath = base_dir + '/AWA014_PTA_2_Rec_Trans/downsampled/'; 
% InputFileTemplate = 'AWA014_PTA_2_';
% MaskPath = ''; 
% ImageOutput = ImagePath + "ST/";
% OutputTemplate = 'AWA014_PTA_2_';
% ImageFileExtension = 'png';
% MaskFileExtension = 'png';
% DataOutput = ImageOutput;
% kindex = 1:256; % updated cropping to make image block as small as possible
% iindex = 1:512;
% jindex = 1:469;
% Nj = length(jindex); Ni = length(iindex); Nk = length(kindex);
% img_digit_size = '%s%s%03d.%s';
% % New variables to use if 
% PermuteOrder = [];
% FlipDim1 = []; % if this is empty then no flip
% FlipDim2 = []; % if this is empty then no flip

%% Full resolution segment version
base_dir = join([getenv("HOME"), "Documents/phd"], '/');
ImagePath = base_dir + '/AWA014_PTA_2_Rec_Trans/'; 
InputFileTemplate = 'AWA014_PTA_2_';
MaskPath = ''; 
ImageOutput = ImagePath + "ST/";
OutputTemplate = 'AWA014_PTA_2_';
ImageFileExtension = 'bmp';
MaskFileExtension = 'png';
DataOutput = ImageOutput;
kindex = 830:1030; % updated cropping to make image block as small as possible
iindex = 220:1220;
jindex = 250:1450;
Nj = length(jindex); Ni = length(iindex); Nk = length(kindex);
img_digit_size = '%s%s%08d.%s';
% New variables to use if 
PermuteOrder = [];
FlipDim1 = []; % if this is empty then no flip
FlipDim2 = []; % if this is empty then no flip

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load mask and image data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(MaskPath)
    fprintf(' ... Loading mask data ... \n');
    M3D = false(Nj,Ni,Nk); 
    parfor k=1:Nk
      fnamein = sprintf(img_digit_size,MaskPath,MaskFileTemplate,kindex(k),MaskFileExtension);
      M = ~(imread(fnamein) == 0);
      M3D(:,:,k) = M(jindex,iindex);
    end
    M3D = permute(M3D,PermuteOrder);
    if ~isempty(FlipDim1)
        M3D = flip(M3D,FlipDim1);
    end
    if ~isempty(FlipDim2)
        M3D = flip(M3D,FlipDim2);
    end
end

fprintf(' ... Loading image data ... \n');
I3D = uint8(false(Nj,Ni,Nk)); 
parfor k=1:Nk
  fnamein = sprintf(img_digit_size,ImagePath,InputFileTemplate,kindex(k),ImageFileExtension);
  I = imread(fnamein);
  I3D(:,:,k) = I(jindex,iindex);
end
%% Commented out the permutation section of the code
% I3D = permute(I3D,PermuteOrder);
% if ~isempty(FlipDim1)
%     I3D = flip(I3D,FlipDim1);
% end
% if ~isempty(FlipDim2)
%     I3D = flip(I3D,FlipDim2);
% end
% 
% % Permute the index arrays also
% ijk = {jindex;iindex;kindex};
% iindex = ijk{PermuteOrder(2)};
% jindex = ijk{PermuteOrder(1)};
% kindex = ijk{PermuteOrder(3)};
% Nj = length(jindex); Ni = length(iindex); Nk = length(kindex);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Mask if necessary
%
% The masking is split into strips that overlap. Each strip is processed
% separately in parallel and then recombined into the final image.
%
% Using this approach the computation is faster as more is being done in
% parallel and it overcomes an issue where a large image stack is not fully
% processed - presumably due to silent memory issues.
%
% Note: split this into two steps/sections - memory issues with one of the
% steps? Do step by step.
%
% Overlap = 100;
% Stride = 600;
% Nk = 1730;
% ProcessStart = max((n-1)*Stride-Overlap+1,1);
% ProcessEnd = min(n*Stride+Overlap,Nk);
% StoreStart = max((n-1)*Stride+1,1);
% StoreEnd = min(n*Stride,Nk);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for strips and working memory allocation
% Overlap = 50;
% Stride = 300;
% Set up for parallel processing. 
Overlap = 10; % after some experimentation, this seems adequate.
% choose stride size so that no more that 12 processes are used. Use 
% ceiling so that stride is at least this number.
NThreads = 12;
Stride = ceil(Nk/NThreads);
Working = cell(ceil(Nk/Stride),1);

% Do processing if is mask is not to be read in
if isempty(MaskPath)
    % Parallel process of independent strips
    parfor i=1:ceil(Nk/Stride)
        Is = max((i-1)*Stride-Overlap+1,1);
        Ie = min(i*Stride+Overlap,Nk);
        fprintf(sprintf('... adjusting and saturating intensity - split %d...\n',i));
        Working{i} = uint8(single(imadjustn(I3D(:,:,Is:Ie))).^1.2);
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
    i = 1;
    StoreEnd = min(Stride,Nk);
    M3D(:,:,1:StoreEnd) = Working{1}(:,:,1:Stride);
    for i=2:ceil(Nk/Stride)-1
        StoreStart = max((i-1)*Stride+1,1);
        StoreEnd = min(i*Stride,Nk);
        M3D(:,:,StoreStart:StoreEnd) = Working{i}(:,:,Overlap+1:end-Overlap);
    end
    StoreStart = max((ceil(Nk/Stride)-1)*Stride+1,1);
    M3D(:,:,StoreStart:Nk) = Working{ceil(Nk/Stride)}(:,:,Overlap+1:end);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Identify ventricle cavity masks by filling cavity and substrating the 
% original mask. This is done as 2D slices in case there is some pathway
% from the cavity to to exterior background in 3D.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('... Iterative 2D filling ...\n');
Working = M3D;
parfor k=1:Nk
    Working(:,:,k) = logical(imfill(Working(:,:,k),'holes'));
end
Working = logical(imsubtract(Working,M3D));

fprintf('... Cleaning ventricular mask ...\n');
stats = regionprops(Working,'Area','PixelIdxList');
[Order,Oloc]=sort(cat(1,stats(:).Area),'descend');
LVMask = false(Nj,Ni,Nk);
LVMask(stats(Oloc(1)).PixelIdxList) = 1; 
RVMask = false(Nj,Ni,Nk);
RVMask(stats(Oloc(2)).PixelIdxList) = 1; 
Mask = LVMask+RVMask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute LV principal component for reslicing images
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find LV lumen mask voxel locations.
fprintf('... Find LV lumen voxel locations ... \n');
stats = regionprops(LVMask,'Area','PixelList','Centroid');
[mv,mloc]=max(cat(1,stats(:).Area));
LVPixelLocs = stats(mloc).PixelList;
LVCentroid = stats(mloc).Centroid;

% Swap the first two coordinates to ensure everything stays in image
% coordiantes - unhelpfully, MATLAB returns the first index of coordinates
% as j and the second index as i. The column swap below corrects this.
LVPixelLocs = LVPixelLocs(:,[2,1,3]);
LVCentroid = LVCentroid(:,[2,1,3]);

fprintf('... Find RV lumen voxel centroid ... \n');
stats = regionprops(RVMask,'Area','Centroid');
[mv,mloc]=max(cat(1,stats(:).Area));
RVCentroid = stats(mloc).Centroid;
RVCentroid = RVCentroid(:,[2,1,3]);

% Use PCA to find main axis of LV voxel cloud
fprintf('... find LV main axis with PCA \n');
Coefs = pca(LVPixelLocs);
a = Coefs(:,1);
% Need to compute -c which is the vector linking the LVCentroid to
% RVCentroid+d*a
% Calculations:
% -c.a = 0 (orthogonal) - note -c as we want the RV centroid to be to the
% left (when viewing the image) of the LV.
% cR + Da + c = cL => c = cL-cR-Da
% therefore: cL.a - cR.a - Da.a = 0
% a.a = 1 since normalised
% D = (cL-cR).a and therefore c = (cL-cR)-((cL-cR).a)a
c = (LVCentroid'-RVCentroid')-(dot((LVCentroid'-RVCentroid'),a))*a;
c = c/norm(c);
b = cross(c,a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code for displaying the ventricles in the material coordinates defined by 
% orthogonal vectors a, b and c.
%
% Note: In resliced image plane b is j and c is i.
%
% Variable XYZConv is required for the reslicing step.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IdxAll = find(M3D);
[XAll,YAll,ZAll] = ind2sub(size(M3D),IdxAll);
T = [a,b,c];
Tdash = inv(T);
XAllDash = XAll-LVCentroid(1);
YAllDash = YAll-LVCentroid(2);
ZAllDash = ZAll-LVCentroid(3);
XYZConv = Tdash*[XAllDash,YAllDash,ZAllDash]';
XYZConv = XYZConv';

IdxLV = find(LVMask);
[XLV,YLV,ZLV] = ind2sub(size(LVMask),IdxLV);
XLVDash = XLV-LVCentroid(1);
YLVDash = YLV-LVCentroid(2);
ZLVDash = ZLV-LVCentroid(3);
XYZLVConv = Tdash*[XLVDash,YLVDash,ZLVDash]';
XYZLVConv = XYZLVConv';

IdxRV = find(RVMask);
[XRV,YRV,ZRV] = ind2sub(size(RVMask),IdxRV);
XRVDash = XRV-LVCentroid(1);
YRVDash = YRV-LVCentroid(2);
ZRVDash = ZRV-LVCentroid(3);
XYZRVConv = Tdash*[XRVDash,YRVDash,ZRVDash]';
XYZRVConv = XYZRVConv';
    
if ~isempty(DataOutput)
    figure;
    scatter3(XYZConv(1:50000:end,1),XYZConv(1:50000:end,2),XYZConv(1:50000:end,3),1,'r','filled'); axis equal; xlabel('{\bfPA}_1'); zlabel('{\bfc}: {\bfRV}_c \rightarrow {\bfLV}_c'); ylabel('{\bfb}={\bfc} \times {\bfPA}_1');
    hold on; scatter3(XYZLVConv(1:50000:end,1),XYZLVConv(1:50000:end,2),XYZLVConv(1:50000:end,3),5,'b','s','filled'); hold off;
    hold on; scatter3(XYZRVConv(1:50000:end,1),XYZRVConv(1:50000:end,2),XYZRVConv(1:50000:end,3),5,'g','s','filled'); hold off;
    view(3); camproj perspective; box on;
    % Print the image for reference
    print(sprintf('%sReorientAxesCoords.png',DataOutput),'-dpng','-r1200');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Find ranges in new image plane, reslice and write images
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ReslicePixelBuffer = 20; % 20 pixel buffer
MinRange = min(round(XYZConv),[],1)-ReslicePixelBuffer; 
MaxRange = max(round(XYZConv),[],1)+ReslicePixelBuffer;
KMin = MinRange(1); KMax = MaxRange(1);
MinRange(1) = 0; MaxRange(1)=0; % The image in the a direction (principal axis) is only 1 voxel thick

if ~isempty(ImageOutput)
    fprintf('...Writing resliced images...\n');
    parfor k=KMin:KMax
        [CImage,Xs,Ys,Zs] = GeneralRotateResample(I3D,LVCentroid+k*a',MinRange,MaxRange,a,b,c);
        fnameout = sprintf('%s%s%03d.png',ImageOutput,OutputTemplate,k-KMin+1);
        imwrite(squeeze(CImage),fnameout);
    end
end


