%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% StructureTensorExtraction.m
%
% This script loads in a stack of images (with isotropic voxels) and 
% computes the components of the structure tensor (the outer product of 
% the intensity gradient vectors). The components are smoothed to a 
% sequence of 1/2 resolutions and these are written to file.
%
% Attempted to save memory by writing intermediate arrays to file and
% by doing calculations in single precision.
%
% Modified by: Mark Trew, October 2014.
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
%ImageClass = 'logical';
%ImageClass = 'uint16'; MaxVal = 2^16;
ImageClass = 'uint8'; MaxVal = 2^8;

% Heart 1 images
Nj = 2048; Ni = 2352; 
kindex = [260:1752];
Nk = length(kindex);
InputFileTemplate = 'FH1_PTA_20_7_21_';
ImageInput = '../ImagesHeart1/MaskedEdited_Extrapolated/';
DataOutput = '../DataHeart1/STBinary/';

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
I = cast(zeros(Nj,Ni,Nk),ImageClass);

% Load in image set
fprintf('... loading images ...\n');
fstring = sprintf('%s%s%%0%dd.%s',ImageInput,InputFileTemplate,DigitsInImageSequence,InputFileExtension);
for k=1:length(kindex),
  if ~mod(k,100), fprintf(' image: %d\n',k); end;
  fnamein = sprintf(fstring,kindex(k)); 
  I(:,:,k) = imread(fnamein);
end;

mu = whos; fprintf(' ** Memory: %f MB\n',sum(cat(1,mu.bytes))/1024/1024);

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
I = permute(I,[2,1,3]);

% Pad the image around each edge - reflect image in padding
fprintf('... padding array ...\n');
IPad = padarray(I,[2,2,2],'symmetric'); clear I;

% Reshape padded permuted equalized image to a 1D array
IPad = reshape(IPad,(Ni+4)*(Nj+4)*(Nk+4),1); 

% convert IPad to a single array between 0 and 1
IPad = single(IPad)/single(max(IPad));

mu = whos; fprintf(' ** Memory: %f MB\n',sum(cat(1,mu.bytes))/1024/1024);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FFT image calculations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute fft of Image
fprintf('... computing FFT of image ...\n');
t10 = clock;
IPadf = fft(IPad); clear IPad;
t1 = etime(clock,t10); fprintf(' fft time for IPad: %0.2f sec\n',t1);
save(strcat(DataOutput,'IPadf.mat'),'-v7.3','IPadf');
clear IPadf;

mu = whos; fprintf(' ** Memory: %f MB\n',sum(cat(1,mu.bytes))/1024/1024);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Derivative filters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the 1D gradient derivative weight filters
fprintf('... computing gradient weight templates ...\n');
[Wx,Wy,Wz] = ConstructDerivativeFilters_SinglePrecision(Ni,Nj,Nk,DerivativeTemplateWidth);
save(strcat(DataOutput,'Wy.mat'),'-v7.3','Wy'); clear Wy;
save(strcat(DataOutput,'Wz.mat'),'-v7.3','Wz'); clear Wz;

mu = whos; fprintf(' ** Memory: %f MB\n',sum(cat(1,mu.bytes))/1024/1024);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FFT derivative calculations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute fft of derivative weights
fprintf('... computing FFTs of derivative weights ...\n');
t10 = clock;
Wxf = fft(Wx); clear Wx; save(strcat(DataOutput,'Wxf.mat'),'-v7.3','Wxf'); clear Wxf;
t1 = etime(clock,t10); fprintf(' FFT time for Wx: %0.2f sec\n',t1);

load(strcat(DataOutput,'Wy.mat'));
t10 = clock;
Wyf = fft(Wy); clear Wy; save(strcat(DataOutput,'Wyf.mat'),'-v7.3','Wyf'); clear Wyf;
t1 = etime(clock,t10); fprintf(' FFT time for Wy: %0.2f sec\n',t1);

load(strcat(DataOutput,'Wz.mat'));
t10 = clock;
Wzf = fft(Wz); clear Wz; % save(strcat(DataOutput,'Wzf.mat'),'-v7.3','Wzf'); % not necessary to save this
t1 = etime(clock,t10); fprintf(' FFT time for Wz: %0.2f sec\n',t1);

% Compute gradients using an fft
fprintf('... performing FFT convolution ...\n');
load(strcat(DataOutput,'IPadf.mat'));
t20 = clock; 
Dk = ifft(IPadf .* Wzf); clear Wzf; save(strcat(DataOutput,'Dk.mat'),'-v7.3','Dk'); clear Dk;% FFT convolution
t2 = etime(clock,t20);
fprintf(' FFT Dz Convolution time: %0.2f sec\n',t2);

load(strcat(DataOutput,'Wxf.mat'));
t20 = clock; 
Di = ifft(IPadf .* Wxf); clear Wxf; save(strcat(DataOutput,'Di.mat'),'-v7.3','Di'); clear Di; % FFT convolution
t2 = etime(clock,t20);
fprintf(' FFT Dx Convolution time: %0.2f sec\n',t2);

load(strcat(DataOutput,'Wyf.mat'));
t20 = clock; 
Dj = ifft(IPadf .* Wyf); clear Wyf; save(strcat(DataOutput,'Dj.mat'),'-v7.3','Dj'); clear Dj;% FFT convolution
t2 = etime(clock,t20);
fprintf(' FFT Dy Convolution time: %0.2f sec\n',t2);

clear IPadf;

% Useful data ranges
PL = 2; % valid derivative padding level
SBi = int16([1+1*PL:Ni+4-1*PL]);
SBj = int16([1+1*PL:Nj+4-1*PL]);
SBk = int16([1+1*PL:Nk+4-1*PL]);

% Raw index locations
[SI0,SJ0,SK0] = ndgrid(int16(1:length(SBi)),int16(1:length(SBj)),int16(1:length(SBk)));

% Report size of zero level data
fprintf('... Zero level data dimensions: (%d,%d,%d)\n',size(SI0));

mu = whos; fprintf(' ** Memory: %f MB\n',sum(cat(1,mu.bytes))/1024/1024);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Perform smoothing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up for first level smoothing
load(strcat(DataOutput,'Di.mat'));
t0 = clock; 
fprintf('... First ST construction and L1 smoothing ...\n');
fprintf('     - Jii ...\n');
Jii = Di.*Di; Jii = reshape(Jii,[(Ni+4),(Nj+4),(Nk+4)]);
fprintf('       class Jii: %s\n',class(Jii));
fprintf('       Level 1\n'); 
[Sii1,SI1,SJ1,SK1] = MultigridAveraging(Jii(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
clear Jii;

fprintf('     - Jij ...\n');
load(strcat(DataOutput,'Dj.mat'));
Jij = Di.*Dj; Jij = reshape(Jij,[(Ni+4),(Nj+4),(Nk+4)]);
[Sij1,SI1,SJ1,SK1] = MultigridAveraging(Jij(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
clear Jij Dj;

fprintf('     - Jik ...\n');
load(strcat(DataOutput,'Dk.mat'));
Jik = Di.*Dk; Jik = reshape(Jik,[(Ni+4),(Nj+4),(Nk+4)]); %clear Di; 
[Sik1,SI1,SJ1,SK1] = MultigridAveraging(Jik(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
clear Jik Di Dk;

fprintf('     - Jjj ...\n');
load(strcat(DataOutput,'Dj.mat'));
Jjj = Dj.*Dj; Jjj = reshape(Jjj,[(Ni+4),(Nj+4),(Nk+4)]);
[Sjj1,SI1,SJ1,SK1] = MultigridAveraging(Jjj(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
clear Jjj;

fprintf('     - Jjk ...\n');
load(strcat(DataOutput,'Dk.mat'));
Jjk = Dj.*Dk; Jjk = reshape(Jjk,[(Ni+4),(Nj+4),(Nk+4)]); %clear Dj; 
[Sjk1,SI1,SJ1,SK1] = MultigridAveraging(Jjk(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
clear Jjk Dj;

fprintf('     - Jkk ...\n');
Jkk = Dk.*Dk; Jkk = reshape(Jkk,[(Ni+4),(Nj+4),(Nk+4)]); %clear Dk; 
[Skk1,SI1,SJ1,SK1] = MultigridAveraging(Jkk(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
clear Jkk;

t1 = etime(clock,t0); fprintf(' First level ST construction and smoothing time: %0.2f sec\n',t1);
fprintf('... First level data dimensions: (%d,%d,%d)\n',size(SI1));

mu = whos; fprintf(' ** Memory: %f MB\n',sum(cat(1,mu.bytes))/1024/1024);

% Smooth to second level using multigrid binomial averaging
fprintf('... Second level smoothing ...\n');
t0 = clock; 
[Sii2,SI2,SJ2,SK2] = MultigridAveraging(Sii1,SI1,SJ1,SK1,SmoothingTemplateWidth);
[Sij2,SI2,SJ2,SK2] = MultigridAveraging(Sij1,SI1,SJ1,SK1,SmoothingTemplateWidth);
[Sik2,SI2,SJ2,SK2] = MultigridAveraging(Sik1,SI1,SJ1,SK1,SmoothingTemplateWidth);
[Sjj2,SI2,SJ2,SK2] = MultigridAveraging(Sjj1,SI1,SJ1,SK1,SmoothingTemplateWidth);
[Sjk2,SI2,SJ2,SK2] = MultigridAveraging(Sjk1,SI1,SJ1,SK1,SmoothingTemplateWidth);
[Skk2,SI2,SJ2,SK2] = MultigridAveraging(Skk1,SI1,SJ1,SK1,SmoothingTemplateWidth);
t1 = etime(clock,t0); fprintf(' Second level smoothing time: %0.2f sec\n',t1);
fprintf('... Second level data dimensions: (%d,%d,%d)\n',size(SI2));

% Smooth to third level using multigrid binomial averaging
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

% Smooth to fourth level using multigrid binomial averaging
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

% Smooth to fifth level using multigrid binomial averaging
fprintf('... Fourth level smoothing ...\n');
t0 = clock; 
[Sii5,SI5,SJ5,SK5] = MultigridAveraging(Sii4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Sij5,SI5,SJ5,SK5] = MultigridAveraging(Sij4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Sik5,SI5,SJ5,SK5] = MultigridAveraging(Sik4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Sjj5,SI5,SJ5,SK5] = MultigridAveraging(Sjj4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Sjk5,SI5,SJ5,SK5] = MultigridAveraging(Sjk4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Skk5,SI5,SJ5,SK5] = MultigridAveraging(Skk4,SI4,SJ4,SK4,SmoothingTemplateWidth);
t1 = etime(clock,t0); fprintf(' Fourth level smoothing time: %0.2f sec\n',t1);
fprintf('... Fifth level data dimensions: (%d,%d,%d)\n',size(SI5));

% Smooth to sixth level using multigrid binomial averaging
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

ttotalprocess1 = cputime; fprintf(' *** total processing time: %0.2f sec\n',ttotalprocess1-ttotalprocess0);

mu = whos; fprintf(' ** Memory: %f MB\n',sum(cat(1,mu.bytes))/1024/1024);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Write binary data files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('... Writing binary data files ...\n');

fid = fopen(sprintf('%sS1.bin',DataOutput),'wb');
fwrite(fid,size(Sii1),'uint16');
fwrite(fid,reshape(Sii1,prod(size(Sii1)),1),'double');
fwrite(fid,reshape(Sij1,prod(size(Sij1)),1),'double');
fwrite(fid,reshape(Sik1,prod(size(Sik1)),1),'double');
fwrite(fid,reshape(Sjj1,prod(size(Sjj1)),1),'double');
fwrite(fid,reshape(Sjk1,prod(size(Sjk1)),1),'double');
fwrite(fid,reshape(Skk1,prod(size(Skk1)),1),'double');
fclose(fid);

fid = fopen(sprintf('%sS2.bin',DataOutput),'wb');
fwrite(fid,size(Sii2),'uint16');
fwrite(fid,reshape(Sii2,prod(size(Sii2)),1),'double');
fwrite(fid,reshape(Sij2,prod(size(Sij2)),1),'double');
fwrite(fid,reshape(Sik2,prod(size(Sik2)),1),'double');
fwrite(fid,reshape(Sjj2,prod(size(Sjj2)),1),'double');
fwrite(fid,reshape(Sjk2,prod(size(Sjk2)),1),'double');
fwrite(fid,reshape(Skk2,prod(size(Skk2)),1),'double');
fclose(fid);

fid = fopen(sprintf('%sS3.bin',DataOutput),'wb');
fwrite(fid,size(Sii3),'uint16');
fwrite(fid,reshape(Sii3,prod(size(Sii3)),1),'double');
fwrite(fid,reshape(Sij3,prod(size(Sij3)),1),'double');
fwrite(fid,reshape(Sik3,prod(size(Sik3)),1),'double');
fwrite(fid,reshape(Sjj3,prod(size(Sjj3)),1),'double');
fwrite(fid,reshape(Sjk3,prod(size(Sjk3)),1),'double');
fwrite(fid,reshape(Skk3,prod(size(Skk3)),1),'double');
fclose(fid);

fid = fopen(sprintf('%sS4.bin',DataOutput),'wb');
fwrite(fid,size(Sii4),'uint16');
fwrite(fid,reshape(Sii4,prod(size(Sii4)),1),'double');
fwrite(fid,reshape(Sij4,prod(size(Sij4)),1),'double');
fwrite(fid,reshape(Sik4,prod(size(Sik4)),1),'double');
fwrite(fid,reshape(Sjj4,prod(size(Sjj4)),1),'double');
fwrite(fid,reshape(Sjk4,prod(size(Sjk4)),1),'double');
fwrite(fid,reshape(Skk4,prod(size(Skk4)),1),'double');
fclose(fid);

fid = fopen(sprintf('%sS5.bin',DataOutput),'wb');
fwrite(fid,size(Sii5),'uint16');
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

fid = fopen(sprintf('%sIJK0.bin',DataOutput),'wb');
fwrite(fid,size(SI0),'uint16');
fwrite(fid,reshape(SI0,prod(size(SI0)),1),'uint16');
fwrite(fid,reshape(SJ0,prod(size(SI0)),1),'uint16');
fwrite(fid,reshape(SK0,prod(size(SI0)),1),'uint16');
fclose(fid);

fid = fopen(sprintf('%sIJK1.bin',DataOutput),'wb');
fwrite(fid,size(SI1),'uint16');
fwrite(fid,reshape(SI1,prod(size(SI1)),1),'uint16');
fwrite(fid,reshape(SJ1,prod(size(SJ1)),1),'uint16');
fwrite(fid,reshape(SK1,prod(size(SK1)),1),'uint16');
fclose(fid);

fid = fopen(sprintf('%sIJK2.bin',DataOutput),'wb');
fwrite(fid,size(SI2),'uint16');
fwrite(fid,reshape(SI2,prod(size(SI2)),1),'uint16');
fwrite(fid,reshape(SJ2,prod(size(SJ2)),1),'uint16');
fwrite(fid,reshape(SK2,prod(size(SK2)),1),'uint16');
fclose(fid);

fid = fopen(sprintf('%sIJK3.bin',DataOutput),'wb');
fwrite(fid,size(SI3),'uint16');
fwrite(fid,reshape(SI3,prod(size(SI3)),1),'uint16');
fwrite(fid,reshape(SJ3,prod(size(SJ3)),1),'uint16');
fwrite(fid,reshape(SK3,prod(size(SK3)),1),'uint16');
fclose(fid);

fid = fopen(sprintf('%sIJK4.bin',DataOutput),'wb');
fwrite(fid,size(SI4),'uint16');
fwrite(fid,reshape(SI4,prod(size(SI4)),1),'uint16');
fwrite(fid,reshape(SJ4,prod(size(SJ4)),1),'uint16');
fwrite(fid,reshape(SK4,prod(size(SK4)),1),'uint16');
fclose(fid);

fid = fopen(sprintf('%sIJK5.bin',DataOutput),'wb');
fwrite(fid,size(SI5),'uint16');
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


