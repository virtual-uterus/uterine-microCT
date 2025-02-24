%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% StructureTensorExtractionTiled.m
%
% This script loads in a stack of images (with isotropic voxels) and 
% computes the components of the structure tensor (the outer product of 
% the intensity gradient vectors). The computations are conducted using
% overlapping tiles spread out memory usage.
% The components are smoothed to a 
% sequence of 1/2 resolutions and these are written to file.
%
% Written by: Mark Trew, October 2021.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

% Timing
ttotalprocess0 = cputime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% User defined parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image class, dimension and numbering sequence
ImageClass = 'uint8'; MaxVal = 2^8;

% % Heart 1 images
% Nj = 2048; Ni = 2352; 
% kindex = [260:1752];
% Nk = length(kindex);
% InputFileTemplate = 'FH1_PTA_20_7_21_';
% ImageInput = '../ImagesHeart1/MaskedEdited_Extrapolated/';
% DataOutput = '../DataHeart1/STBinary/';

% Heart Air Dry images
% Nj = 1237; Ni = 1237;
% kindex = [548:1060];
% Nk = length(kindex);
% InputFileTemplate = 'AirDry_';
% ImageInput = '../ImagesHeartDry/MaskedExtrapolated/';
% DataOutput = '../DataHeartDry/STBinary/';

% H1C1H images
% Nj = 1589; Ni = 1771;
% kindex = [370:1689];
% Nk = length(kindex);
% InputFileTemplate = 'H1C1H_';
% ImageInput = '../Images_H1C1H/MaskedCleanExtrapolated/';
% DataOutput = '../Data_H1C1H/STBinary/';

% H1FGR1H images
% Nj = 1586; Ni = 2083;
% kindex = [1:1460];
% Nk = length(kindex);
% InputFileTemplate = 'H1FGR1H_';
% ImageInput = '../Images_H1FGR1H/MaskedCleanExtrapolated/';
% DataOutput = '../Images_H1FGR1H/Data_H1FGR1H/STBinary/';

% H1C1H images 20 Sept 2022
Nj = 1613; Ni = 2118;
kindex = [1:1616];
Nk = length(kindex);
InputFileTemplate = 'H1C1H_';
ImageInput = '../Images_H1C1H/MaskedCleanExtrapolated/';
DataOutput = '../Images_H1C1H/Data_H1C1H/STBinary/';

% Other parameteres
InputFileExtension = 'png';
DigitsInImageSequence = 4; % number of digits in image numbering pattern

% Set the derivative and smoothing template voxel widths
DerivativeTemplateWidth = 5;
%DerivativeTemplateWidth = 3;
SmoothingTemplateWidth = 5;
%SmoothingTemplateWidth = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set up and load in image stack
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory
fprintf('... allocating memory ...\n');
ID = cast(zeros(Nj,Ni,Nk),ImageClass);

% Load in image set
fprintf('... loading images ...\n');
fstring = sprintf('%s%s%%0%dd.%s',ImageInput,InputFileTemplate,DigitsInImageSequence,InputFileExtension);
for k=1:length(kindex)
  if ~mod(k,100), fprintf(' image: %d\n',k); end
  fnamein = sprintf(fstring,kindex(k)); 
  ID(:,:,k) = imread(fnamein);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Rearrange and pad the image stack
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Permute the array to make the x/i index the most rapidly varying.
% This process is necessary to ensure the i,j,k indices of the matlab
% array correspond with the xyz coordinates of the texture orientation
% vectors.
fprintf('... permuting array ...\n');
ID = permute(ID,[2,1,3]);
ID = single(ID)./single(max(ID(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Break up image and iterate over sections
%
% NSi - number of whole blocks of width Si in image
% Si - width of NSi blocks
% SiStar - width of remaining partial block
% 
% Block widths are chosen to make the partial block width as close to the 
% other block widths as possible.
% 
% The image is padded with a width of W on all sides.
%
% The actual image subblocks are Si+2W etc in size.
%
% The position of the starting corner of the subblock in absolute indices
% of padded array is (for subblock (I,J,K)):
%   [(I-1)*Si+1,(J-1)*Sj,(K-1)*Sk]
% The position of the finishing corner of the subblock is:
%   [I*Si+2*W,J*Sj+2*W,K*Sk+2*W] for I=1:NSi, J=1:NSj, K=1:NSk
%   [I*SiStar+2*W,J*SjStar+2*W,K*SkStar+2*W] for I=NSi+1, J=NSj+1, K=NSk+1
%
% The relative positions of the unpadded image subblocks (that tesselate
% to give the original image):
%   Start: [W+1,W+1,W+1]
%   End:   [Si+W,Sj+W,Sk+W] for I=1:NSi, J=1:NSj, K=1:NSk
%          [SiStar+W,SjStar+W,SkStar+W] for I=NSi+1, J=NSj+1, K=NSk+1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('... striding array ...\n');

W = 2; % the template width is 5 for DTW = 5 or 3
OL = 2*W; % Overlap is wider than template to ensure all points are supported

% Pad array with overlap to make extraction of useful information easier
ID = padarray(ID,[W,W,W],'symmetric','both');
[Nip,Njp,Nkp]=size(ID);

% Compute a good stride
S = [100:1:200]; % limits on S %300
NSi = floor(Ni./S); Pi = S-(Ni-NSi.*S);
NSj = floor(Nj./S); Pj = S-(Nj-NSj.*S);
NSk = floor(Nk./S); Pk = S-(Nk-NSk.*S);
[MOi,IMOi]=min((Pi/max(Pi))+(NSi/max(NSi))); NSi = NSi(IMOi); Si = S(IMOi); SiStar = Ni-NSi*Si;
[MOj,IMOj]=min((Pj/max(Pj))+(NSj/max(NSj))); NSj = NSj(IMOj); Sj = S(IMOj); SjStar = Nj-NSj*Sj;
[MOk,IMOk]=min((Pk/max(Pk))+(NSk/max(NSk))); NSk = NSk(IMOk); Sk = S(IMOk); SkStar = Nk-NSk*Sk;

% % Write overlap information
% fprintf('Array size: (%d,%d,%d), Overlap: %d\n',Ni,Nj,Nk,OL);
fprintf('Number of tiles: (%d,%d,%d)\n',NSi,NSj,NSk);
fprintf('Tile widths: (%d,%d,%d)\n',Si,Sj,Sk);
% 
% Loop over subunits
t10 = clock;

% Level 4 binomial filter
Wm = [1,4,6,4,1]/16; 
PL = 2; % padding level for filter

% Double step ranges and offset for smoothing
I = [1+1*PL:2:Ni-1*PL]; Nmi = length(I); MaxI = max(I)+2;
J = [1+1*PL:2:Nj-1*PL]; Nmj = length(J); MaxJ = max(J)+2;

% Storage of a row in j
ST = cell(6,1);
for n=1:6
    ST{n} = cast(false(Nmi+2,Nmj+2,NSj*Sk+SkStar),class(ID));
end
    
% Repeated k units
for K=1:NSk
    fprintf('K: %d\n',K);
    Temp = TiledJUnits(ID,K,NSi,NSj,[Si,Sj,Sk],[Si,Sj,Sk],SiStar,SjStar,W,DerivativeTemplateWidth);    
    for n=1:6
         ST{n}(:,:,(K-1)*Sk+1:(K-1)*Sk+Sk) = Temp{n};
    end
end
% Last k unit with k width SkStar
Temp = TiledJUnits(ID,NSk+1,NSi,NSj,[Si,Sj,Sk],[Si,Sj,SkStar],SiStar,SjStar,W,DerivativeTemplateWidth);
for n=1:6
    ST{n}(:,:,NSk*Sk+1:NSk*Sk+SkStar,:) = Temp{n};
end
clear Temp;

% Smooth in k
fprintf(' ... Smoothing in k.');
K = [1+1*PL:2:Nk-1*PL]; Nmk = length(K); MaxK = max(K)+2;
K1 = [3,K-2,MaxK-2];
K2 = [2,K-1,MaxK-1];
K3 = [1,K,MaxK];
K4 = [2,K+1,MaxK-1];
K5 = [3,K+2,MaxK-2];
for n=1:6
    ST{n} = Wm(1)*ST{n}(:,:,K1)+Wm(2)*ST{n}(:,:,K2)+Wm(3)*ST{n}(:,:,K3)+Wm(4)*ST{n}(:,:,K4)+Wm(5)*ST{n}(:,:,K5);
end
fprintf(' size 2: (%d,%d,%d)\n',size(ST{1}));    

t1 = etime(clock,t10); fprintf(' Total tile time: %0.2f sec\n',t1);
 
% Perform one level of extraction on indices (L1) and then do all smoothing
% together (L2-L7) and then write binary files
% ST is already one level smoothed
% Raw index locations
SI1 = single([1,I,MaxI]);
SJ1 = single([1,J,MaxJ]);
SK1 = single([1,K,MaxK]);
[SI1,SJ1,SK1] = ndgrid(SI1,SJ1,SK1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smooth to second level using multigrid binomial averaging
fprintf('... Second level smoothing ...\n');
t0 = clock; 
[Sii2,SI2,SJ2,SK2] = MultigridAveraging(ST{1},SI1,SJ1,SK1,SmoothingTemplateWidth);
[Sij2,SI2,SJ2,SK2] = MultigridAveraging(ST{2},SI1,SJ1,SK1,SmoothingTemplateWidth);
[Sik2,SI2,SJ2,SK2] = MultigridAveraging(ST{3},SI1,SJ1,SK1,SmoothingTemplateWidth);
[Sjj2,SI2,SJ2,SK2] = MultigridAveraging(ST{4},SI1,SJ1,SK1,SmoothingTemplateWidth);
[Sjk2,SI2,SJ2,SK2] = MultigridAveraging(ST{5},SI1,SJ1,SK1,SmoothingTemplateWidth);
[Skk2,SI2,SJ2,SK2] = MultigridAveraging(ST{6},SI1,SJ1,SK1,SmoothingTemplateWidth);
t1 = etime(clock,t0); fprintf(' Second level smoothing time: %0.2f sec\n',t1);
fprintf('... Second level data dimensions: (%d,%d,%d)\n',size(SI2));

fprintf('... Third level smoothing ...\n');
t0 = clock; 
[Sii3,SI3,SJ3,SK3] = MultigridAveraging(Sii2,SI2,SJ2,SK2,SmoothingTemplateWidth);
[Sij3,SI3,SJ3,SK3] = MultigridAveraging(Sij2,SI2,SJ2,SK2,SmoothingTemplateWidth);
[Sik3,SI3,SJ3,SK3] = MultigridAveraging(Sik2,SI2,SJ2,SK2,SmoothingTemplateWidth);
[Sjj3,SI3,SJ3,SK3] = MultigridAveraging(Sjj2,SI2,SJ2,SK2,SmoothingTemplateWidth);
[Sjk3,SI3,SJ3,SK3] = MultigridAveraging(Sjk2,SI2,SJ2,SK2,SmoothingTemplateWidth);
[Skk3,SI3,SJ3,SK3] = MultigridAveraging(Skk2,SI2,SJ2,SK2,SmoothingTemplateWidth);
t1 = etime(clock,t0); fprintf(' Third level smoothing time: %0.2f sec\n',t1);
fprintf('... Third level data dimensions: (%d,%d,%d)\n',size(SI3));

fprintf('... Fourth level smoothing ...\n');
t0 = clock; 
[Sii4,SI4,SJ4,SK4] = MultigridAveraging(Sii3,SI3,SJ3,SK3,SmoothingTemplateWidth);
[Sij4,SI4,SJ4,SK4] = MultigridAveraging(Sij3,SI3,SJ3,SK3,SmoothingTemplateWidth);
[Sik4,SI4,SJ4,SK4] = MultigridAveraging(Sik3,SI3,SJ3,SK3,SmoothingTemplateWidth);
[Sjj4,SI4,SJ4,SK4] = MultigridAveraging(Sjj3,SI3,SJ3,SK3,SmoothingTemplateWidth);
[Sjk4,SI4,SJ4,SK4] = MultigridAveraging(Sjk3,SI3,SJ3,SK3,SmoothingTemplateWidth);
[Skk4,SI4,SJ4,SK4] = MultigridAveraging(Skk3,SI3,SJ3,SK3,SmoothingTemplateWidth);
t1 = etime(clock,t0); fprintf(' Fourth level smoothing time: %0.2f sec\n',t1);
fprintf('... Fourth level data dimensions: (%d,%d,%d)\n',size(SI4));

fprintf('... Fifth level smoothing ...\n');
t0 = clock; 
[Sii5,SI5,SJ5,SK5] = MultigridAveraging(Sii4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Sij5,SI5,SJ5,SK5] = MultigridAveraging(Sij4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Sik5,SI5,SJ5,SK5] = MultigridAveraging(Sik4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Sjj5,SI5,SJ5,SK5] = MultigridAveraging(Sjj4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Sjk5,SI5,SJ5,SK5] = MultigridAveraging(Sjk4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Skk5,SI5,SJ5,SK5] = MultigridAveraging(Skk4,SI4,SJ4,SK4,SmoothingTemplateWidth);
t1 = etime(clock,t0); fprintf(' Fifth level smoothing time: %0.2f sec\n',t1);
fprintf('... Fifth level data dimensions: (%d,%d,%d)\n',size(SI5));

fprintf('... Sixth level smoothing ...\n');
t0 = clock; 
[Sii6,SI6,SJ6,SK6] = MultigridAveraging(Sii5,SI5,SJ5,SK5,SmoothingTemplateWidth);
[Sij6,SI6,SJ6,SK6] = MultigridAveraging(Sij5,SI5,SJ5,SK5,SmoothingTemplateWidth);
[Sik6,SI6,SJ6,SK6] = MultigridAveraging(Sik5,SI5,SJ5,SK5,SmoothingTemplateWidth);
[Sjj6,SI6,SJ6,SK6] = MultigridAveraging(Sjj5,SI5,SJ5,SK5,SmoothingTemplateWidth);
[Sjk6,SI6,SJ6,SK6] = MultigridAveraging(Sjk5,SI5,SJ5,SK5,SmoothingTemplateWidth);
[Skk6,SI6,SJ6,SK6] = MultigridAveraging(Skk5,SI5,SJ5,SK5,SmoothingTemplateWidth);
t1 = etime(clock,t0); fprintf(' Sixth level smoothing time: %0.2f sec\n',t1);
fprintf('... Sixth level data dimensions: (%d,%d,%d)\n',size(SI6));

fprintf('... Seventh level smoothing ...\n');
t0 = clock; 
[Sii7,SI7,SJ7,SK7] = MultigridAveraging(Sii6,SI6,SJ6,SK6,SmoothingTemplateWidth);
[Sij7,SI7,SJ7,SK7] = MultigridAveraging(Sij6,SI6,SJ6,SK6,SmoothingTemplateWidth);
[Sik7,SI7,SJ7,SK7] = MultigridAveraging(Sik6,SI6,SJ6,SK6,SmoothingTemplateWidth);
[Sjj7,SI7,SJ7,SK7] = MultigridAveraging(Sjj6,SI6,SJ6,SK6,SmoothingTemplateWidth);
[Sjk7,SI7,SJ7,SK7] = MultigridAveraging(Sjk6,SI6,SJ6,SK6,SmoothingTemplateWidth);
[Skk7,SI7,SJ7,SK7] = MultigridAveraging(Skk6,SI6,SJ6,SK6,SmoothingTemplateWidth);
t1 = etime(clock,t0); fprintf(' Seventh level smoothing time: %0.2f sec\n',t1);
fprintf('... Seventh level data dimensions: (%d,%d,%d)\n',size(SI7));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing binary data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('... Writing binary data files ...\n');

fid = fopen(sprintf('%sS2.bin',DataOutput),'wb');
fwrite(fid,[size(Sii2,1),size(Sii2,2),size(Sii2,3)],'uint16');
fwrite(fid,reshape(Sii2,prod(size(Sii2)),1),'double');
fwrite(fid,reshape(Sij2,prod(size(Sij2)),1),'double');
fwrite(fid,reshape(Sik2,prod(size(Sik2)),1),'double');
fwrite(fid,reshape(Sjj2,prod(size(Sjj2)),1),'double');
fwrite(fid,reshape(Sjk2,prod(size(Sjk2)),1),'double');
fwrite(fid,reshape(Skk2,prod(size(Skk2)),1),'double');
fclose(fid);

fid = fopen(sprintf('%sS3.bin',DataOutput),'wb');
fwrite(fid,[size(Sii3,1),size(Sii3,2),size(Sii3,3)],'uint16');
fwrite(fid,reshape(Sii3,prod(size(Sii3)),1),'double');
fwrite(fid,reshape(Sij3,prod(size(Sij3)),1),'double');
fwrite(fid,reshape(Sik3,prod(size(Sik3)),1),'double');
fwrite(fid,reshape(Sjj3,prod(size(Sjj3)),1),'double');
fwrite(fid,reshape(Sjk3,prod(size(Sjk3)),1),'double');
fwrite(fid,reshape(Skk3,prod(size(Skk3)),1),'double');
fclose(fid);

fid = fopen(sprintf('%sS4.bin',DataOutput),'wb');
fwrite(fid,[size(Sii4,1),size(Sii4,2),size(Sii4,3)],'uint16');
fwrite(fid,reshape(Sii4,prod(size(Sii4)),1),'double');
fwrite(fid,reshape(Sij4,prod(size(Sij4)),1),'double');
fwrite(fid,reshape(Sik4,prod(size(Sik4)),1),'double');
fwrite(fid,reshape(Sjj4,prod(size(Sjj4)),1),'double');
fwrite(fid,reshape(Sjk4,prod(size(Sjk4)),1),'double');
fwrite(fid,reshape(Skk4,prod(size(Skk4)),1),'double');
fclose(fid);

fid = fopen(sprintf('%sS5.bin',DataOutput),'wb');
fwrite(fid,[size(Sii5,1),size(Sii5,2),size(Sii5,3)],'uint16');
fwrite(fid,reshape(Sii5,prod(size(Sii5)),1),'double');
fwrite(fid,reshape(Sij5,prod(size(Sij5)),1),'double');
fwrite(fid,reshape(Sik5,prod(size(Sik5)),1),'double');
fwrite(fid,reshape(Sjj5,prod(size(Sjj5)),1),'double');
fwrite(fid,reshape(Sjk5,prod(size(Sjk5)),1),'double');
fwrite(fid,reshape(Skk5,prod(size(Skk5)),1),'double');
fclose(fid);

fid = fopen(sprintf('%sS6.bin',DataOutput),'wb');
fwrite(fid,[size(Sii6,1),size(Sii6,2),size(Sii6,3)],'uint16');
fwrite(fid,reshape(Sii6,prod(size(Sii6)),1),'double');
fwrite(fid,reshape(Sij6,prod(size(Sij6)),1),'double');
fwrite(fid,reshape(Sik6,prod(size(Sik6)),1),'double');
fwrite(fid,reshape(Sjj6,prod(size(Sjj6)),1),'double');
fwrite(fid,reshape(Sjk6,prod(size(Sjk6)),1),'double');
fwrite(fid,reshape(Skk6,prod(size(Skk6)),1),'double');
fclose(fid);

fid = fopen(sprintf('%sS7.bin',DataOutput),'wb');
fwrite(fid,[size(Sii7,1),size(Sii7,2),size(Sii7,3)],'uint16');
fwrite(fid,reshape(Sii7,prod(size(Sii7)),1),'double');
fwrite(fid,reshape(Sij7,prod(size(Sij7)),1),'double');
fwrite(fid,reshape(Sik7,prod(size(Sik7)),1),'double');
fwrite(fid,reshape(Sjj7,prod(size(Sjj7)),1),'double');
fwrite(fid,reshape(Sjk7,prod(size(Sjk7)),1),'double');
fwrite(fid,reshape(Skk7,prod(size(Skk7)),1),'double');
fclose(fid);

fid = fopen(sprintf('%sIJK2.bin',DataOutput),'wb');
fwrite(fid,[size(SI2,1),size(SI2,2),size(SI2,3)],'uint16');
fwrite(fid,reshape(SI2,prod(size(SI2)),1),'uint16');
fwrite(fid,reshape(SJ2,prod(size(SJ2)),1),'uint16');
fwrite(fid,reshape(SK2,prod(size(SK2)),1),'uint16');
fclose(fid);

fid = fopen(sprintf('%sIJK3.bin',DataOutput),'wb');
fwrite(fid,[size(SI3,1),size(SI3,2),size(SI3,3)],'uint16');
fwrite(fid,reshape(SI3,prod(size(SI3)),1),'uint16');
fwrite(fid,reshape(SJ3,prod(size(SJ3)),1),'uint16');
fwrite(fid,reshape(SK3,prod(size(SK3)),1),'uint16');
fclose(fid);

fid = fopen(sprintf('%sIJK4.bin',DataOutput),'wb');
fwrite(fid,[size(SI4,1),size(SI4,2),size(SI4,3)],'uint16');
fwrite(fid,reshape(SI4,prod(size(SI4)),1),'uint16');
fwrite(fid,reshape(SJ4,prod(size(SJ4)),1),'uint16');
fwrite(fid,reshape(SK4,prod(size(SK4)),1),'uint16');
fclose(fid);

fid = fopen(sprintf('%sIJK5.bin',DataOutput),'wb');
fwrite(fid,[size(SI5,1),size(SI5,2),size(SI5,3)],'uint16');
fwrite(fid,reshape(SI5,prod(size(SI5)),1),'uint16');
fwrite(fid,reshape(SJ5,prod(size(SJ5)),1),'uint16');
fwrite(fid,reshape(SK5,prod(size(SK5)),1),'uint16');
fclose(fid);

fid = fopen(sprintf('%sIJK6.bin',DataOutput),'wb');
fwrite(fid,[size(SI6,1),size(SI6,2),size(SI6,3)],'uint16');
fwrite(fid,reshape(SI6,prod(size(SI6)),1),'uint16');
fwrite(fid,reshape(SJ6,prod(size(SJ6)),1),'uint16');
fwrite(fid,reshape(SK6,prod(size(SK6)),1),'uint16');
fclose(fid);

fid = fopen(sprintf('%sIJK7.bin',DataOutput),'wb');
fwrite(fid,[size(SI7,1),size(SI7,2),size(SI7,3)],'uint16');
fwrite(fid,reshape(SI7,prod(size(SI7)),1),'uint16');
fwrite(fid,reshape(SJ7,prod(size(SJ7)),1),'uint16');
fwrite(fid,reshape(SK7,prod(size(SK7)),1),'uint16');
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ST = TiledJUnits(ID,K,NSi,NSj,SS,SE,SiStar,SjStar,W,DTW)

    % Level 4 binomial filter
    Wm = [1,4,6,4,1]/16; 
    PL = 2; % padding level for filter

    % Double step ranges and offset for smoothing
    I = [1+1*PL:2:(NSi*SS(1)+SiStar)-1*PL]; Nmi = length(I);
    J = [1+1*PL:2:(NSj*SS(2)+SjStar)-1*PL]; Nmj = length(J); MaxJ = max(J)+2;

    % Storage of a row in j
    ST = cell(6,1);
    for n=1:6
        ST{n} = cast(false(Nmi+2,NSj*SS(2)+SjStar,SE(3)),class(ID));
    end
    
    % Repeated j units
    for J=1:NSj        
        fprintf('J: %d\n',J);
        Temp = TiledIUnits(ID,J,K,NSi,SS,SE,SiStar,W,DTW);
        for n=1:6
             ST{n}(:,(J-1)*SS(2)+1:(J-1)*SS(2)+SE(2),:) = Temp{n};
        end
    end
    Temp = TiledIUnits(ID,NSj+1,K,NSi,SS,[SE(1),SjStar,SE(3)],SiStar,W,DTW);
    for n=1:6
        ST{n}(:,NSj*SS(2)+1:NSj*SS(2)+SjStar,:) = Temp{n};
    end
    
    % Smooth in j
    fprintf(' ... Smoothing in j.');
    fprintf(' size 1: (%d,%d,%d)',size(ST{1}));
    J = [1+1*PL:2:(NSj*SS(2)+SjStar)-1*PL]; Nmj = length(J); MaxJ = max(J)+2;
    J1 = [3,J-2,MaxJ-2];
    J2 = [2,J-1,MaxJ-1];
    J3 = [1,J,MaxJ];
    J4 = [2,J+1,MaxJ-1];
    J5 = [3,J+2,MaxJ-2];
    for n=1:6
        ST{n} = Wm(1)*ST{n}(:,J1,:)+Wm(2)*ST{n}(:,J2,:)+Wm(3)*ST{n}(:,J3,:)+Wm(4)*ST{n}(:,J4,:)+Wm(5)*ST{n}(:,J5,:);
    end
    fprintf(' size 2: (%d,%d,%d)\n',size(ST{1}));    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ST = TiledIUnits(ID,J,K,NSi,SS,SE,SiStar,W,DTW)

    % This function tiles units in I
    % 
    % SS - width of unit for start calc
    % SE - width of unit for end calc (same as SS except for last one)
    % 

    % Full storage of a row in i
    ST = cell(6,1);
    for n=1:6
        ST{n} = cast(false(NSi*SS(1)+SiStar,SE(2),SE(3)),class(ID));
    end
    
    % Repeated i units
    for I=1:NSi
        IS = ID((I-1)*SS(1)+1:(I-1)*SS(1)+SE(1)+2*W,(J-1)*SS(2)+1:(J-1)*SS(2)+SE(2)+2*W,(K-1)*SS(3)+1:(K-1)*SS(3)+SE(3)+2*W);
        Temp = ComputeStructureTensor(IS,DTW);
        for n=1:6
            ST{n}((I-1)*SS(1)+1:(I-1)*SS(1)+SE(1),:,:) = Temp{n}(W+1:W+SE(1),W+1:W+SE(2),W+1:W+SE(3));
        end
    end
    % Last i unit
    IS = ID((NSi)*SS(1)+1:NSi*SS(1)+SiStar+2*W,(J-1)*SS(2)+1:(J-1)*SS(2)+SE(2)+2*W,(K-1)*SS(3)+1:(K-1)*SS(3)+SE(3)+2*W);
    Temp = ComputeStructureTensor(IS,DTW);
    for n=1:6
        ST{n}(NSi*SS(1)+1:NSi*SS(1)+SiStar,:,:) = Temp{n}(W+1:W+SiStar,W+1:W+SE(2),W+1:W+SE(3));
    end
    
    % Level 4 binomial filter
    Wm = [1,4,6,4,1]/16; 
    PL = 2; % padding level for filter

    % Double step ranges and offset for smoothing
    I = [1+1*PL:2:(NSi*SS(1)+SiStar)-1*PL]; Nmi = length(I); MaxI = max(I)+2;
  
    % Smooth in i
    fprintf(' ... Smoothing in i.');
    fprintf(' size 1: (%d,%d,%d)',size(ST{1}));
    I1 = [3,I-2,MaxI-2];
    I2 = [2,I-1,MaxI-1];
    I3 = [1,I,MaxI];
    I4 = [2,I+1,MaxI-1];
    I5 = [3,I+2,MaxI-2];
    for n=1:6
        ST{n} = Wm(1)*ST{n}(I1,:,:)+Wm(2)*ST{n}(I2,:,:)+Wm(3)*ST{n}(I3,:,:)+Wm(4)*ST{n}(I4,:,:)+Wm(5)*ST{n}(I5,:,:);
    end
    fprintf(' size 2: (%d,%d,%d)\n',size(ST{1}));   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

