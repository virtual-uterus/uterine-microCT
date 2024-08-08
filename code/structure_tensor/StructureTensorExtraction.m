%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% StructureTensorExtraction.m
%
% This script loads in a stack of images (with isotropic voxels) and 
% computes the components of the structure tensor (the outer product of 
% the intensity gradient vectors). The components are smoothed to a 
% sequence of 1/2 resolutions and these are written to file.
%
% Modified by: Mark Trew, September 2021.
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

%%
kindex = 80:534; 
Nj = 548; Ni = 1124;
Nk = length(kindex);
base_dir = join([getenv("HOME"), "Documents/phd"], '/');
ImageInput = base_dir + '/AWA015_PTA_1_Rec_Trans/downsampled/ST/extrapolated/';
InputFileTemplate = 'AWA015_PTA_1_';
DataOutput = base_dir + '/AWA015_PTA_1_Rec_Trans/downsampled/ST/binary/';
img_digit_size = '%s%s%03d.%s';
%%
% kindex = 1:325; 
% Nj = 518; Ni = 753;
% Nk = length(kindex);
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% ImageInput = base_dir + '/AWA015_PTA_2_Ova_Rec_Trans/downsampled/ST/extrapolated/';
% InputFileTemplate = 'AWA015_PTA_2_Ova_Trans_';
% DataOutput = base_dir + '/AWA015_PTA_2_Ova_Rec_Trans/downsampled/ST/binary/';
% img_digit_size = '%s%s%03d.%s';
%% Downsampled version
% Nj = 469; Ni = 512;
% kindex = 1:138;
% Nk = length(kindex);
% InputFileTemplate = 'AWA014_PTA_2_';
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% ImageInput = base_dir + '/AWA014_PTA_2_Rec_Trans/downsampled/ST/extrapolated/';
% DataOutput = base_dir + '/AWA014_PTA_2_Rec_Trans/downsampled/ST/binary/';

%% Full resolution segment version
% kindex = 830:1030; 
% Nj = 1401; Ni = 1001;
% Nk = length(kindex);
% InputFileTemplate = 'AWA014_PTA_2_';
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% ImageInput = base_dir + '/AWA014_PTA_2_Rec_Trans/ST/extrapolated/';
% DataOutput = base_dir + '/AWA014_PTA_2_Rec_Trans/ST/binary/';

%% Reoriented downsampled version
% Nj = 367; Ni = 473;
% kindex = 1:352;
% Nk = length(kindex);
% InputFileTemplate = 'AWA014_PTA_2_';
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% ImageInput = base_dir + '/AWA014_PTA_2_Rec_Trans/downsampled/ST/extrapolated/';
% DataOutput = base_dir + '/AWA014_PTA_2_Rec_Trans/downsampled/ST/binary/';

%% Reoriented full resolution segment version
% Nj = 420; Ni = 388;
% kindex = 1:383;
% Nk = length(kindex);
% InputFileTemplate = 'AWA014_PTA_2_';
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% ImageInput = base_dir + '/AWA014_PTA_2_Rec_Trans/ST/extrapolated/';
% DataOutput = base_dir + '/AWA014_PTA_2_Rec_Trans/ST/binary/';


% Other parameteres
InputFileExtension = 'png';
DigitsInImageSequence = 3; % number of digits in image numbering pattern

% Set the derivative and smoothing template voxel widths
DerivativeTemplateWidth = 3;
%DerivativeTemplateWidth = 3;
SmoothingTemplateWidth = 3;
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
for k=1:length(kindex)
  if ~mod(k,100), fprintf(' image: %d\n',k); end
  fnamein = sprintf(fstring,kindex(k)); 
  I(:,:,k) = imread(fnamein);
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
I = permute(I,[2,1,3]);

% Pad the image around each edge - reflect image in padding
fprintf('... padding array ...\n');
IPad = padarray(I,[2,2,2],'symmetric'); clear I;

% Reshape padded permuted equalized image to a 1D array
% IPadO = IPad;
IPad = reshape(IPad,(Ni+4)*(Nj+4)*(Nk+4),1); 

% convert IPad to a double array between 0 and 1
% IPadO = double(IPadO)./double(max(IPadO(:)));
%IPad = double(IPad)/double(max(IPad));
IPad = single(IPad)/single(max(IPad));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Derivative filters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the 1D gradient derivative weight filters
fprintf('... computing gradient weight templates ...\n');
[Wx,Wy,Wz] = ConstructDerivativeFilters(Ni,Nj,Nk,DerivativeTemplateWidth);

% Set up the 1D 2nd derivative weight filters
% fprintf('... computing 2nd derivative weight templates ...\n');
% [Wxx,Wyy,Wzz,Wxy,Wxz,Wyz] = Construct2ndDerivativeFilters(Ni,Nj,Nk,DerivativeTemplateWidth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FFT calculations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute fft of Image
fprintf('... computing FFT of image ...\n');
t10 = clock;
IPadf = single(fft(IPad)); %clear IPad;
t1 = etime(clock,t10); fprintf(' fft time for IPad: %0.2f sec\n',t1);

% % Compute fft of 1st derivative weights
% fprintf('... computing FFTs of 1st derivative weights ...\n');
% t10 = clock;
% Wxf = fft(Wx); clear Wx;
% t1 = etime(clock,t10); fprintf(' FFT time for Wx: %0.2f sec\n',t1);
% t10 = clock;
% Wyf = fft(Wy); clear Wy;
% t1 = etime(clock,t10); fprintf(' FFT time for Wy: %0.2f sec\n',t1);
% t10 = clock;
% Wzf = fft(Wz); clear Wz;
% t1 = etime(clock,t10); fprintf(' FFT time for Wz: %0.2f sec\n',t1);

% Compute fft of 2nd derivative weights
% fprintf('... computing FFTs of 2nd derivative weights ...\n');
% t10 = clock;
% Wxxf = fft(Wxx); clear Wxx;
% t1 = etime(clock,t10); fprintf(' FFT time for Wxx: %0.2f sec\n',t1);
% t10 = clock;
% Wyyf = fft(Wyy); clear Wyy;
% t1 = etime(clock,t10); fprintf(' FFT time for Wyy: %0.2f sec\n',t1);
% t10 = clock;
% Wzzf = fft(Wzz); clear Wzz;
% t1 = etime(clock,t10); fprintf(' FFT time for Wzz: %0.2f sec\n',t1);
% t10 = clock;
% Wxyf = fft(Wxy); clear Wxy;
% t1 = etime(clock,t10); fprintf(' FFT time for Wxy: %0.2f sec\n',t1);
% t10 = clock;
% Wxzf = fft(Wxz); clear Wxz;
% t1 = etime(clock,t10); fprintf(' FFT time for Wxz: %0.2f sec\n',t1);
% t10 = clock;
% Wyzf = fft(Wyz); clear Wyz;
% t1 = etime(clock,t10); fprintf(' FFT time for Wyz: %0.2f sec\n',t1);

% Compute fft of 1st derivative weights
fprintf('... computing FFTs of 1st derivative weights ...\n');
fprintf('... performing FFT convolution for 1st derivatives ...\n');
t10 = clock;
Wxf = single(fft(Wx)); clear Wx;
Di = double(ifft(IPadf .* Wxf)); clear Wxf;% FFT convolution
t1 = etime(clock,t10); fprintf(' FFT Dx Convolution time for Wx: %0.2f sec\n',t1);
t10 = clock;
Wyf = single(fft(Wy)); clear Wy;
Dj = double(ifft(IPadf .* Wyf)); clear Wyf;% FFT convolution
t1 = etime(clock,t10); fprintf(' FFT Dy Convolution time for Wy: %0.2f sec\n',t1);
t10 = clock;
Wzf = single(fft(Wz)); clear Wz;
Dk = double(ifft(IPadf .* Wzf)); clear Wzf;% FFT convolution
t1 = etime(clock,t10); fprintf(' FFT Dz Convolution time for Wz: %0.2f sec\n',t1);

% Compute gradients using an fft
% fprintf('... performing FFT convolution for 1st derivatives ...\n');
% t20 = clock; 
% Di = ifft(IPadf .* Wxf); clear Wxf;% FFT convolution
% t2 = etime(clock,t20);
% fprintf(' FFT Dx Convolution time: %0.2f sec\n',t2);
% t20 = clock; 
% Dj = ifft(IPadf .* Wyf); clear Wyf;% FFT convolution
% t2 = etime(clock,t20);
% fprintf(' FFT Dy Convolution time: %0.2f sec\n',t2);
% t20 = clock; 
% Dk = ifft(IPadf .* Wzf); clear Wzf;% FFT convolution
% t2 = etime(clock,t20);
% fprintf(' FFT Dz Convolution time: %0.2f sec\n',t2);

% Std gradient calcs
% fprintf('... Using inbuilt gradient function ...\n');
% [Dj,Di,Dk] = gradient(IPadO);
% Di = reshape(Di,[(Ni+4)*(Nj+4)*(Nk+4),1]);
% Dj = reshape(Dj,[(Ni+4)*(Nj+4)*(Nk+4),1]);
% Dk = reshape(Dk,[(Ni+4)*(Nj+4)*(Nk+4),1]);

% Compute second derivatives using an fft
% fprintf('... performing FFT convolution for 2nd derivatives ...\n');
% t20 = clock; 
% Hii = ifft(IPadf .* Wxxf); clear Wxxf;% FFT convolution
% t2 = etime(clock,t20);
% fprintf(' FFT Hxx Convolution time: %0.2f sec\n',t2);
% t20 = clock; 
% Hjj = ifft(IPadf .* Wyyf); clear Wyyf;% FFT convolution
% t2 = etime(clock,t20);
% fprintf(' FFT Hyy Convolution time: %0.2f sec\n',t2);
% t20 = clock; 
% Hkk = ifft(IPadf .* Wzzf); clear Wzzf;% FFT convolution
% t2 = etime(clock,t20);
% fprintf(' FFT Hzz Convolution time: %0.2f sec\n',t2);
% t20 = clock; 
% Hij = ifft(IPadf .* Wxyf); clear Wxyf;% FFT convolution
% t2 = etime(clock,t20);
% fprintf(' FFT Hxy Convolution time: %0.2f sec\n',t2);
% t20 = clock; 
% Hik = ifft(IPadf .* Wxzf); clear Wxzf;% FFT convolution
% t2 = etime(clock,t20);
% fprintf(' FFT Hxz Convolution time: %0.2f sec\n',t2);
% t20 = clock; 
% Hjk = ifft(IPadf .* Wyzf); clear Wyzf;% FFT convolution
% t2 = etime(clock,t20);
% fprintf(' FFT Hyz Convolution time: %0.2f sec\n',t2);

clear IPadf;

% Useful data ranges
PL = 2; % valid derivative padding level
SBi = [1+1*PL:Ni+4-1*PL];
SBj = [1+1*PL:Nj+4-1*PL];
SBk = [1+1*PL:Nk+4-1*PL];

% Raw index locations
[SI0,SJ0,SK0] = ndgrid(1:length(SBi),1:length(SBj),1:length(SBk));

% Report size of zero level data
fprintf('... Zero level data dimensions: (%d,%d,%d)\n',size(SI0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Perform smoothing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up for first level smoothing
t0 = clock; 
fprintf('... First level ST construction and smoothing ...\n');

% fprintf('     - Dii ...\n');
% Dii = Di; Dii = reshape(Dii,[(Ni+4),(Nj+4),(Nk+4)]); 
% [Dii1,SI1,SJ1,SK1] = MultigridAveraging(Dii(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
% 
% fprintf('     - Djj ...\n');
% Djj = Dj; Djj = reshape(Djj,[(Ni+4),(Nj+4),(Nk+4)]); 
% [Djj1,SI1,SJ1,SK1] = MultigridAveraging(Djj(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
% 
% fprintf('     - Dkk ...\n');
% Dkk = Dk; Dkk = reshape(Dkk,[(Ni+4),(Nj+4),(Nk+4)]); 
% [Dkk1,SI1,SJ1,SK1] = MultigridAveraging(Dkk(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
% 
fprintf('     - Jii ...\n');
Jii = Di.*Di; Jii = reshape(Jii,[(Ni+4),(Nj+4),(Nk+4)]); 
[Sii1,SI1,SJ1,SK1] = MultigridAveraging(Jii(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);

fprintf('     - Jij ...\n');
Jij = Di.*Dj; Jij = reshape(Jij,[(Ni+4),(Nj+4),(Nk+4)]);
[Sij1,SI1,SJ1,SK1] = MultigridAveraging(Jij(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);

fprintf('     - Jik ...\n');
Jik = Di.*Dk; Jik = reshape(Jik,[(Ni+4),(Nj+4),(Nk+4)]); %clear Di; 
[Sik1,SI1,SJ1,SK1] = MultigridAveraging(Jik(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);

fprintf('     - Jjj ...\n');
Jjj = Dj.*Dj; Jjj = reshape(Jjj,[(Ni+4),(Nj+4),(Nk+4)]);
[Sjj1,SI1,SJ1,SK1] = MultigridAveraging(Jjj(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);

fprintf('     - Jjk ...\n');
Jjk = Dj.*Dk; Jjk = reshape(Jjk,[(Ni+4),(Nj+4),(Nk+4)]); %clear Dj; 
[Sjk1,SI1,SJ1,SK1] = MultigridAveraging(Jjk(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);

fprintf('     - Jkk ...\n');
Jkk = Dk.*Dk; Jkk = reshape(Jkk,[(Ni+4),(Nj+4),(Nk+4)]); %clear Dk; 
[Skk1,SI1,SJ1,SK1] = MultigridAveraging(Jkk(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);

t1 = etime(clock,t0); fprintf(' First level ST construction and smoothing time: %0.2f sec\n',t1);

% t0 = clock; 
% fprintf('... First level H construction and smoothing ...\n');
% fprintf('     - Hii ...\n');
% Hii = reshape(Hii,[(Ni+4),(Nj+4),(Nk+4)]); 
% [Hii1,SI1,SJ1,SK1] = MultigridAveraging(Hii(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
% fprintf('     - Hjj ...\n');
% Hjj = reshape(Hjj,[(Ni+4),(Nj+4),(Nk+4)]); 
% [Hjj1,SI1,SJ1,SK1] = MultigridAveraging(Hjj(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
% fprintf('     - Hkk ...\n');
% Hkk = reshape(Hkk,[(Ni+4),(Nj+4),(Nk+4)]); 
% [Hkk1,SI1,SJ1,SK1] = MultigridAveraging(Hkk(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
% fprintf('     - Hij ...\n');
% Hij = reshape(Hij,[(Ni+4),(Nj+4),(Nk+4)]);
% [Hij1,SI1,SJ1,SK1] = MultigridAveraging(Hij(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
% fprintf('     - Hik ...\n');
% Hik = reshape(Hik,[(Ni+4),(Nj+4),(Nk+4)]); %clear Di; 
% [Hik1,SI1,SJ1,SK1] = MultigridAveraging(Hik(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
% fprintf('     - Hjk ...\n');
% Hjk = reshape(Hjk,[(Ni+4),(Nj+4),(Nk+4)]); %clear Dj; 
% [Hjk1,SI1,SJ1,SK1] = MultigridAveraging(Hjk(SBi,SBj,SBk),SI0,SJ0,SK0,SmoothingTemplateWidth);
% t1 = etime(clock,t0); fprintf(' First level H construction and smoothing time: %0.2f sec\n',t1);

fprintf('... First level data dimensions: (%d,%d,%d)\n',size(SI1));

% Smooth to second level using multigrid binomial averaging
fprintf('... Second level smoothing ...\n');
t0 = clock; 
% [Dii2,SI2,SJ2,SK2] = MultigridAveraging(Dii1,SI1,SJ1,SK1,SmoothingTemplateWidth);
% [Djj2,SI2,SJ2,SK2] = MultigridAveraging(Djj1,SI1,SJ1,SK1,SmoothingTemplateWidth);
% [Dkk2,SI2,SJ2,SK2] = MultigridAveraging(Dkk1,SI1,SJ1,SK1,SmoothingTemplateWidth);
[Sii2,SI2,SJ2,SK2] = MultigridAveraging(Sii1,SI1,SJ1,SK1,SmoothingTemplateWidth);
[Sij2,SI2,SJ2,SK2] = MultigridAveraging(Sij1,SI1,SJ1,SK1,SmoothingTemplateWidth);
[Sik2,SI2,SJ2,SK2] = MultigridAveraging(Sik1,SI1,SJ1,SK1,SmoothingTemplateWidth);
[Sjj2,SI2,SJ2,SK2] = MultigridAveraging(Sjj1,SI1,SJ1,SK1,SmoothingTemplateWidth);
[Sjk2,SI2,SJ2,SK2] = MultigridAveraging(Sjk1,SI1,SJ1,SK1,SmoothingTemplateWidth);
[Skk2,SI2,SJ2,SK2] = MultigridAveraging(Skk1,SI1,SJ1,SK1,SmoothingTemplateWidth);
% [Hii2,SI2,SJ2,SK2] = MultigridAveraging(Hii1,SI1,SJ1,SK1,SmoothingTemplateWidth);
% [Hij2,SI2,SJ2,SK2] = MultigridAveraging(Hij1,SI1,SJ1,SK1,SmoothingTemplateWidth);
% [Hik2,SI2,SJ2,SK2] = MultigridAveraging(Hik1,SI1,SJ1,SK1,SmoothingTemplateWidth);
% [Hjj2,SI2,SJ2,SK2] = MultigridAveraging(Hjj1,SI1,SJ1,SK1,SmoothingTemplateWidth);
% [Hjk2,SI2,SJ2,SK2] = MultigridAveraging(Hjk1,SI1,SJ1,SK1,SmoothingTemplateWidth);
% [Hkk2,SI2,SJ2,SK2] = MultigridAveraging(Hkk1,SI1,SJ1,SK1,SmoothingTemplateWidth);
t1 = etime(clock,t0); fprintf(' Second level smoothing time: %0.2f sec\n',t1);
fprintf('... Second level data dimensions: (%d,%d,%d)\n',size(SI2));

% Smooth to third level using multigrid binomial averaging
fprintf('... Third level smoothing ...\n');
t0 = clock; 
% [Dii3,SI3,SJ3,SK3] = MultigridAveraging(Dii2,SI2,SJ2,SK2,SmoothingTemplateWidth);
% [Djj3,SI3,SJ3,SK3] = MultigridAveraging(Djj2,SI2,SJ2,SK2,SmoothingTemplateWidth);
% [Dkk3,SI3,SJ3,SK3] = MultigridAveraging(Dkk2,SI2,SJ2,SK2,SmoothingTemplateWidth);
[Sii3,SI3,SJ3,SK3] = MultigridAveraging(Sii2,SI2,SJ2,SK2,SmoothingTemplateWidth);
[Sij3,SI3,SJ3,SK3] = MultigridAveraging(Sij2,SI2,SJ2,SK2,SmoothingTemplateWidth);
[Sik3,SI3,SJ3,SK3] = MultigridAveraging(Sik2,SI2,SJ2,SK2,SmoothingTemplateWidth);
[Sjj3,SI3,SJ3,SK3] = MultigridAveraging(Sjj2,SI2,SJ2,SK2,SmoothingTemplateWidth);
[Sjk3,SI3,SJ3,SK3] = MultigridAveraging(Sjk2,SI2,SJ2,SK2,SmoothingTemplateWidth);
[Skk3,SI3,SJ3,SK3] = MultigridAveraging(Skk2,SI2,SJ2,SK2,SmoothingTemplateWidth);
% [Hii3,SI3,SJ3,SK3] = MultigridAveraging(Hii2,SI2,SJ2,SK2,SmoothingTemplateWidth);
% [Hij3,SI3,SJ3,SK3] = MultigridAveraging(Hij2,SI2,SJ2,SK2,SmoothingTemplateWidth);
% [Hik3,SI3,SJ3,SK3] = MultigridAveraging(Hik2,SI2,SJ2,SK2,SmoothingTemplateWidth);
% [Hjj3,SI3,SJ3,SK3] = MultigridAveraging(Hjj2,SI2,SJ2,SK2,SmoothingTemplateWidth);
% [Hjk3,SI3,SJ3,SK3] = MultigridAveraging(Hjk2,SI2,SJ2,SK2,SmoothingTemplateWidth);
% [Hkk3,SI3,SJ3,SK3] = MultigridAveraging(Hkk2,SI2,SJ2,SK2,SmoothingTemplateWidth);
t1 = etime(clock,t0); fprintf(' Third level smoothing time: %0.2f sec\n',t1);
fprintf('... Third level data dimensions: (%d,%d,%d)\n',size(SI3));

% Smooth to fourth level using multigrid binomial averaging
fprintf('... Fourth level smoothing ...\n');
t0 = clock; 
% [Dii4,SI4,SJ4,SK4] = MultigridAveraging(Dii3,SI3,SJ3,SK3,SmoothingTemplateWidth);
% [Djj4,SI4,SJ4,SK4] = MultigridAveraging(Djj3,SI3,SJ3,SK3,SmoothingTemplateWidth);
% [Dkk4,SI4,SJ4,SK4] = MultigridAveraging(Dkk3,SI3,SJ3,SK3,SmoothingTemplateWidth);
[Sii4,SI4,SJ4,SK4] = MultigridAveraging(Sii3,SI3,SJ3,SK3,SmoothingTemplateWidth);
[Sij4,SI4,SJ4,SK4] = MultigridAveraging(Sij3,SI3,SJ3,SK3,SmoothingTemplateWidth);
[Sik4,SI4,SJ4,SK4] = MultigridAveraging(Sik3,SI3,SJ3,SK3,SmoothingTemplateWidth);
[Sjj4,SI4,SJ4,SK4] = MultigridAveraging(Sjj3,SI3,SJ3,SK3,SmoothingTemplateWidth);
[Sjk4,SI4,SJ4,SK4] = MultigridAveraging(Sjk3,SI3,SJ3,SK3,SmoothingTemplateWidth);
[Skk4,SI4,SJ4,SK4] = MultigridAveraging(Skk3,SI3,SJ3,SK3,SmoothingTemplateWidth);
% [Hii4,SI4,SJ4,SK4] = MultigridAveraging(Hii3,SI3,SJ3,SK3,SmoothingTemplateWidth);
% [Hij4,SI4,SJ4,SK4] = MultigridAveraging(Hij3,SI3,SJ3,SK3,SmoothingTemplateWidth);
% [Hik4,SI4,SJ4,SK4] = MultigridAveraging(Hik3,SI3,SJ3,SK3,SmoothingTemplateWidth);
% [Hjj4,SI4,SJ4,SK4] = MultigridAveraging(Hjj3,SI3,SJ3,SK3,SmoothingTemplateWidth);
% [Hjk4,SI4,SJ4,SK4] = MultigridAveraging(Hjk3,SI3,SJ3,SK3,SmoothingTemplateWidth);
% [Hkk4,SI4,SJ4,SK4] = MultigridAveraging(Hkk3,SI3,SJ3,SK3,SmoothingTemplateWidth);
t1 = etime(clock,t0); fprintf(' Fourth level smoothing time: %0.2f sec\n',t1);
fprintf('... Fourth level data dimensions: (%d,%d,%d)\n',size(SI4));

% Smooth to fifth level using multigrid binomial averaging
fprintf('... Fifth level smoothing ...\n');
t0 = clock; 
% [Dii5,SI5,SJ5,SK5] = MultigridAveraging(Dii4,SI4,SJ4,SK4,SmoothingTemplateWidth);
% [Djj5,SI5,SJ5,SK5] = MultigridAveraging(Djj4,SI4,SJ4,SK4,SmoothingTemplateWidth);
% [Dkk5,SI5,SJ5,SK5] = MultigridAveraging(Dkk4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Sii5,SI5,SJ5,SK5] = MultigridAveraging(Sii4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Sij5,SI5,SJ5,SK5] = MultigridAveraging(Sij4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Sik5,SI5,SJ5,SK5] = MultigridAveraging(Sik4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Sjj5,SI5,SJ5,SK5] = MultigridAveraging(Sjj4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Sjk5,SI5,SJ5,SK5] = MultigridAveraging(Sjk4,SI4,SJ4,SK4,SmoothingTemplateWidth);
[Skk5,SI5,SJ5,SK5] = MultigridAveraging(Skk4,SI4,SJ4,SK4,SmoothingTemplateWidth);
% [Hii5,SI5,SJ5,SK5] = MultigridAveraging(Hii4,SI4,SJ4,SK4,SmoothingTemplateWidth);
% [Hij5,SI5,SJ5,SK5] = MultigridAveraging(Hij4,SI4,SJ4,SK4,SmoothingTemplateWidth);
% [Hik5,SI5,SJ5,SK5] = MultigridAveraging(Hik4,SI4,SJ4,SK4,SmoothingTemplateWidth);
% [Hjj5,SI5,SJ5,SK5] = MultigridAveraging(Hjj4,SI4,SJ4,SK4,SmoothingTemplateWidth);
% [Hjk5,SI5,SJ5,SK5] = MultigridAveraging(Hjk4,SI4,SJ4,SK4,SmoothingTemplateWidth);
% [Hkk5,SI5,SJ5,SK5] = MultigridAveraging(Hkk4,SI4,SJ4,SK4,SmoothingTemplateWidth);
t1 = etime(clock,t0); fprintf(' Fifth level smoothing time: %0.2f sec\n',t1);
fprintf('... Fifth level data dimensions: (%d,%d,%d)\n',size(SI5));

% Smooth to sixth level using multigrid binomial averaging
fprintf('... Sixth level smoothing ...\n');
t0 = clock; 
% [Dii6,SI6,SJ6,SK6] = MultigridAveraging(Dii5,SI5,SJ5,SK5,SmoothingTemplateWidth);
% [Djj6,SI6,SJ6,SK6] = MultigridAveraging(Djj5,SI5,SJ5,SK5,SmoothingTemplateWidth);
% [Dkk6,SI6,SJ6,SK6] = MultigridAveraging(Dkk5,SI5,SJ5,SK5,SmoothingTemplateWidth);
[Sii6,SI6,SJ6,SK6] = MultigridAveraging(Sii5,SI5,SJ5,SK5,SmoothingTemplateWidth);
[Sij6,SI6,SJ6,SK6] = MultigridAveraging(Sij5,SI5,SJ5,SK5,SmoothingTemplateWidth);
[Sik6,SI6,SJ6,SK6] = MultigridAveraging(Sik5,SI5,SJ5,SK5,SmoothingTemplateWidth);
[Sjj6,SI6,SJ6,SK6] = MultigridAveraging(Sjj5,SI5,SJ5,SK5,SmoothingTemplateWidth);
[Sjk6,SI6,SJ6,SK6] = MultigridAveraging(Sjk5,SI5,SJ5,SK5,SmoothingTemplateWidth);
[Skk6,SI6,SJ6,SK6] = MultigridAveraging(Skk5,SI5,SJ5,SK5,SmoothingTemplateWidth);
% [Hii6,SI6,SJ6,SK6] = MultigridAveraging(Hii5,SI5,SJ5,SK5,SmoothingTemplateWidth);
% [Hij6,SI6,SJ6,SK6] = MultigridAveraging(Hij5,SI5,SJ5,SK5,SmoothingTemplateWidth);
% [Hik6,SI6,SJ6,SK6] = MultigridAveraging(Hik5,SI5,SJ5,SK5,SmoothingTemplateWidth);
% [Hjj6,SI6,SJ6,SK6] = MultigridAveraging(Hjj5,SI5,SJ5,SK5,SmoothingTemplateWidth);
% [Hjk6,SI6,SJ6,SK6] = MultigridAveraging(Hjk5,SI5,SJ5,SK5,SmoothingTemplateWidth);
% [Hkk6,SI6,SJ6,SK6] = MultigridAveraging(Hkk5,SI5,SJ5,SK5,SmoothingTemplateWidth);

% Smooth to seventh level using multigrid binomial averaging
fprintf('... Seventh level smoothing ...\n');
t0 = clock; 
% [Dii7,SI7,SJ7,SK7] = MultigridAveraging(Dii6,SI6,SJ6,SK6,SmoothingTemplateWidth);
% [Djj7,SI7,SJ7,SK7] = MultigridAveraging(Djj6,SI6,SJ6,SK6,SmoothingTemplateWidth);
% [Dkk7,SI7,SJ7,SK7] = MultigridAveraging(Dkk6,SI6,SJ6,SK6,SmoothingTemplateWidth);
[Sii7,SI7,SJ7,SK7] = MultigridAveraging(Sii6,SI6,SJ6,SK6,SmoothingTemplateWidth);
[Sij7,SI7,SJ7,SK7] = MultigridAveraging(Sij6,SI6,SJ6,SK6,SmoothingTemplateWidth);
[Sik7,SI7,SJ7,SK7] = MultigridAveraging(Sik6,SI6,SJ6,SK6,SmoothingTemplateWidth);
[Sjj7,SI7,SJ7,SK7] = MultigridAveraging(Sjj6,SI6,SJ6,SK6,SmoothingTemplateWidth);
[Sjk7,SI7,SJ7,SK7] = MultigridAveraging(Sjk6,SI6,SJ6,SK6,SmoothingTemplateWidth);
[Skk7,SI7,SJ7,SK7] = MultigridAveraging(Skk6,SI6,SJ6,SK6,SmoothingTemplateWidth);
% [Hii7,SI7,SJ7,SK7] = MultigridAveraging(Hii6,SI6,SJ6,SK6,SmoothingTemplateWidth);
% [Hij7,SI7,SJ7,SK7] = MultigridAveraging(Hij6,SI6,SJ6,SK6,SmoothingTemplateWidth);
% [Hik7,SI7,SJ7,SK7] = MultigridAveraging(Hik6,SI6,SJ6,SK6,SmoothingTemplateWidth);
% [Hjj7,SI7,SJ7,SK7] = MultigridAveraging(Hjj6,SI6,SJ6,SK6,SmoothingTemplateWidth);
% [Hjk7,SI7,SJ7,SK7] = MultigridAveraging(Hjk6,SI6,SJ6,SK6,SmoothingTemplateWidth);
% [Hkk7,SI7,SJ7,SK7] = MultigridAveraging(Hkk6,SI6,SJ6,SK6,SmoothingTemplateWidth);
t1 = etime(clock,t0); fprintf(' Seventh level smoothing time: %0.2f sec\n',t1);
fprintf('... Seventh level data dimensions: (%d,%d,%d)\n',size(SI7));

ttotalprocess1 = cputime; fprintf(' *** total processing time: %0.2f sec\n',ttotalprocess1-ttotalprocess0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Write binary data files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('... Writing binary data files ...\n');

% fid = fopen(sprintf('%sS0.bin',DataOutput),'wb');
% fwrite(fid,[size(Jii,1),size(Jii,2),size(Jii,3)],'uint16');
% fwrite(fid,reshape(Jii,prod(size(Jii)),1),'double');
% fwrite(fid,reshape(Jij,prod(size(Jij)),1),'double');
% fwrite(fid,reshape(Jik,prod(size(Jik)),1),'double');
% fwrite(fid,reshape(Jjj,prod(size(Jjj)),1),'double');
% fwrite(fid,reshape(Jjk,prod(size(Jjk)),1),'double');
% fwrite(fid,reshape(Jkk,prod(size(Jkk)),1),'double');
% % fwrite(fid,reshape(Dii,prod(size(Dii)),1),'double');
% % fwrite(fid,reshape(Djj,prod(size(Djj)),1),'double');
% % fwrite(fid,reshape(Dkk,prod(size(Dkk)),1),'double');
% fclose(fid);

% fid = fopen(sprintf('%sS1.bin',DataOutput),'wb');
% fwrite(fid,[size(Sii1,1),size(Sii1,2),size(Sii1,3)],'uint16');
% fwrite(fid,reshape(Sii1,prod(size(Sii1)),1),'double');
% fwrite(fid,reshape(Sij1,prod(size(Sij1)),1),'double');
% fwrite(fid,reshape(Sik1,prod(size(Sik1)),1),'double');
% fwrite(fid,reshape(Sjj1,prod(size(Sjj1)),1),'double');
% fwrite(fid,reshape(Sjk1,prod(size(Sjk1)),1),'double');
% fwrite(fid,reshape(Skk1,prod(size(Skk1)),1),'double');
% % fwrite(fid,reshape(Dii1,prod(size(Sii1)),1),'double');
% % fwrite(fid,reshape(Djj1,prod(size(Sii1)),1),'double');
% % fwrite(fid,reshape(Dkk1,prod(size(Sii1)),1),'double');
% fclose(fid);

fid = fopen(sprintf('%sS2.bin',DataOutput),'wb');
fwrite(fid,[size(Sii2,1),size(Sii2,2),size(Sii2,3)],'uint16');
fwrite(fid,reshape(Sii2,numel(Sii2),1),'double');
fwrite(fid,reshape(Sij2,numel(Sij2),1),'double');
fwrite(fid,reshape(Sik2,numel(Sik2),1),'double');
fwrite(fid,reshape(Sjj2,numel(Sjj2),1),'double');
fwrite(fid,reshape(Sjk2,numel(Sjk2),1),'double');
fwrite(fid,reshape(Skk2,numel(Skk2),1),'double');
% fwrite(fid,reshape(Dii2,prod(size(Sii2)),1),'double');
% fwrite(fid,reshape(Djj2,prod(size(Sii2)),1),'double');
% fwrite(fid,reshape(Dkk2,prod(size(Sii2)),1),'double');
fclose(fid);

fid = fopen(sprintf('%sS3.bin',DataOutput),'wb');
fwrite(fid,[size(Sii3,1),size(Sii3,2),size(Sii3,3)],'uint16');
fwrite(fid,reshape(Sii3,numel(Sii3),1),'double');
fwrite(fid,reshape(Sij3,numel(Sij3),1),'double');
fwrite(fid,reshape(Sik3,numel(Sik3),1),'double');
fwrite(fid,reshape(Sjj3,numel(Sjj3),1),'double');
fwrite(fid,reshape(Sjk3,numel(Sjk3),1),'double');
fwrite(fid,reshape(Skk3,numel(Skk3),1),'double');
% fwrite(fid,reshape(Dii3,prod(size(Sii3)),1),'double');
% fwrite(fid,reshape(Djj3,prod(size(Sii3)),1),'double');
% fwrite(fid,reshape(Dkk3,prod(size(Sii3)),1),'double');
fclose(fid);

fid = fopen(sprintf('%sS4.bin',DataOutput),'wb');
fwrite(fid,[size(Sii4,1),size(Sii4,2),size(Sii4,3)],'uint16');
fwrite(fid,reshape(Sii4,numel(Sii4),1),'double');
fwrite(fid,reshape(Sij4,numel(Sij4),1),'double');
fwrite(fid,reshape(Sik4,numel(Sik4),1),'double');
fwrite(fid,reshape(Sjj4,numel(Sjj4),1),'double');
fwrite(fid,reshape(Sjk4,numel(Sjk4),1),'double');
fwrite(fid,reshape(Skk4,numel(Skk4),1),'double');
% fwrite(fid,reshape(Dii4,prod(size(Sii4)),1),'double');
% fwrite(fid,reshape(Djj4,prod(size(Sii4)),1),'double');
% fwrite(fid,reshape(Dkk4,prod(size(Sii4)),1),'double');
fclose(fid);

fid = fopen(sprintf('%sS5.bin',DataOutput),'wb');
fwrite(fid,[size(Sii5,1),size(Sii5,2),size(Sii5,3)],'uint16');
fwrite(fid,reshape(Sii5,numel(Sii5),1),'double');
fwrite(fid,reshape(Sij5,numel(Sij5),1),'double');
fwrite(fid,reshape(Sik5,numel(Sik5),1),'double');
fwrite(fid,reshape(Sjj5,numel(Sjj5),1),'double');
fwrite(fid,reshape(Sjk5,numel(Sjk5),1),'double');
fwrite(fid,reshape(Skk5,numel(Skk5),1),'double');
% fwrite(fid,reshape(Dii5,prod(size(Sii5)),1),'double');
% fwrite(fid,reshape(Djj5,prod(size(Sii5)),1),'double');
% fwrite(fid,reshape(Dkk5,prod(size(Sii5)),1),'double');
% fwrite(fid,reshape(Hii5,prod(size(Sii5)),1),'double');
% fwrite(fid,reshape(Hij5,prod(size(Sij5)),1),'double');
% fwrite(fid,reshape(Hik5,prod(size(Sik5)),1),'double');
% fwrite(fid,reshape(Hjj5,prod(size(Sjj5)),1),'double');
% fwrite(fid,reshape(Hjk5,prod(size(Sjk5)),1),'double');
% fwrite(fid,reshape(Hkk5,prod(size(Skk5)),1),'double');
fclose(fid);

fid = fopen(sprintf('%sS6.bin',DataOutput),'wb');
fwrite(fid,[size(Sii6,1),size(Sii6,2),size(Sii6,3)],'uint16');
fwrite(fid,reshape(Sii6,numel(Sii6),1),'double');
fwrite(fid,reshape(Sij6,numel(Sij6),1),'double');
fwrite(fid,reshape(Sik6,numel(Sik6),1),'double');
fwrite(fid,reshape(Sjj6,numel(Sjj6),1),'double');
fwrite(fid,reshape(Sjk6,numel(Sjk6),1),'double');
fwrite(fid,reshape(Skk6,numel(Skk6),1),'double');
% fwrite(fid,reshape(Dii6,prod(size(Sii6)),1),'double');
% fwrite(fid,reshape(Djj6,prod(size(Sii6)),1),'double');
% fwrite(fid,reshape(Dkk6,prod(size(Sii6)),1),'double');
% fwrite(fid,reshape(Hii6,prod(size(Sii6)),1),'double');
% fwrite(fid,reshape(Hij6,prod(size(Sij6)),1),'double');
% fwrite(fid,reshape(Hik6,prod(size(Sik6)),1),'double');
% fwrite(fid,reshape(Hjj6,prod(size(Sjj6)),1),'double');
% fwrite(fid,reshape(Hjk6,prod(size(Sjk6)),1),'double');
% fwrite(fid,reshape(Hkk6,prod(size(Skk6)),1),'double');
fclose(fid);

fid = fopen(sprintf('%sS7.bin',DataOutput),'wb');
fwrite(fid,[size(Sii7,1),size(Sii7,2),size(Sii7,3)],'uint16');
fwrite(fid,reshape(Sii7,numel(Sii7),1),'double');
fwrite(fid,reshape(Sij7,numel(Sij7),1),'double');
fwrite(fid,reshape(Sik7,numel(Sik7),1),'double');
fwrite(fid,reshape(Sjj7,numel(Sjj7),1),'double');
fwrite(fid,reshape(Sjk7,numel(Sjk7),1),'double');
fwrite(fid,reshape(Skk7,numel(Skk7),1),'double');
% fwrite(fid,reshape(Dii7,prod(size(Sii7)),1),'double');
% fwrite(fid,reshape(Djj7,prod(size(Sii7)),1),'double');
% fwrite(fid,reshape(Dkk7,prod(size(Sii7)),1),'double');
% fwrite(fid,reshape(Hii7,prod(size(Sii7)),1),'double');
% fwrite(fid,reshape(Hij7,prod(size(Sij7)),1),'double');
% fwrite(fid,reshape(Hik7,prod(size(Sik7)),1),'double');
% fwrite(fid,reshape(Hjj7,prod(size(Sjj7)),1),'double');
% fwrite(fid,reshape(Hjk7,prod(size(Sjk7)),1),'double');
% fwrite(fid,reshape(Hkk7,prod(size(Skk7)),1),'double');
fclose(fid);

% fid = fopen(sprintf('%sIJK0.bin',DataOutput),'wb');
% fwrite(fid,[size(SI0,1),size(SI0,2),size(SI0,3)],'uint16');
% fwrite(fid,reshape(SI0,prod(size(SI0)),1),'uint16');
% fwrite(fid,reshape(SJ0,prod(size(SI0)),1),'uint16');
% fwrite(fid,reshape(SK0,prod(size(SI0)),1),'uint16');
% fclose(fid);
% 
% fid = fopen(sprintf('%sIJK1.bin',DataOutput),'wb');
% fwrite(fid,[size(SI1,1),size(SI1,2),size(SI1,3)],'uint16');
% fwrite(fid,reshape(SI1,prod(size(SI1)),1),'uint16');
% fwrite(fid,reshape(SJ1,prod(size(SJ1)),1),'uint16');
% fwrite(fid,reshape(SK1,prod(size(SK1)),1),'uint16');
% fclose(fid);

fid = fopen(sprintf('%sIJK2.bin',DataOutput),'wb');
fwrite(fid,[size(SI2,1),size(SI2,2),size(SI2,3)],'uint16');
fwrite(fid,reshape(SI2,numel(SI2),1),'uint16');
fwrite(fid,reshape(SJ2,numel(SJ2),1),'uint16');
fwrite(fid,reshape(SK2,numel(SK2),1),'uint16');
fclose(fid);

fid = fopen(sprintf('%sIJK3.bin',DataOutput),'wb');
fwrite(fid,[size(SI3,1),size(SI3,2),size(SI3,3)],'uint16');
fwrite(fid,reshape(SI3,numel(SI3),1),'uint16');
fwrite(fid,reshape(SJ3,numel(SJ3),1),'uint16');
fwrite(fid,reshape(SK3,numel(SK3),1),'uint16');
fclose(fid);

fid = fopen(sprintf('%sIJK4.bin',DataOutput),'wb');
fwrite(fid,[size(SI4,1),size(SI4,2),size(SI4,3)],'uint16');
fwrite(fid,reshape(SI4,numel(SI4),1),'uint16');
fwrite(fid,reshape(SJ4,numel(SJ4),1),'uint16');
fwrite(fid,reshape(SK4,numel(SK4),1),'uint16');
fclose(fid);

fid = fopen(sprintf('%sIJK5.bin',DataOutput),'wb');
fwrite(fid,[size(SI5,1),size(SI5,2),size(SI5,3)],'uint16');
fwrite(fid,reshape(SI5,numel(SI5),1),'uint16');
fwrite(fid,reshape(SJ5,numel(SJ5),1),'uint16');
fwrite(fid,reshape(SK5,numel(SK5),1),'uint16');
fclose(fid);

fid = fopen(sprintf('%sIJK6.bin',DataOutput),'wb');
fwrite(fid,[size(SI6,1),size(SI6,2),size(SI6,3)],'uint16');
fwrite(fid,reshape(SI6,numel(SI6),1),'uint16');
fwrite(fid,reshape(SJ6,numel(SJ6),1),'uint16');
fwrite(fid,reshape(SK6,numel(SK6),1),'uint16');
fclose(fid);

fid = fopen(sprintf('%sIJK7.bin',DataOutput),'wb');
fwrite(fid,[size(SI7,1),size(SI7,2),size(SI7,3)],'uint16');
fwrite(fid,reshape(SI7,numel(SI7),1),'uint16');
fwrite(fid,reshape(SJ7,numel(SJ7),1),'uint16');
fwrite(fid,reshape(SK7,numel(SK7),1),'uint16');
fclose(fid);



