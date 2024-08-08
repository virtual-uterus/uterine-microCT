%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ProcessST_MaskHelix.m
%
% This script loads in structure tensor binary data at a specific level
% of smoothing, tissue mask images and information about the centroid position and
% long axis of the ventricles.
%
% Modified by: Mark Trew, June 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

ttotalprocess0 = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set parameters 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Level = 5; % frequency resolution of ST/Hessian data to use

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Data and path locations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full image set - heart 2 original
% InputPath = '../Data/STBinary/';
% OutputPath = '../Data/STBinary/';
% MaskPath = '../RawImages/BinaryMask/';
% HelixImageOutputPath = '../Images/Part1HelixAngleImages/Level5/';
% MaskPath = '../Images/Mask/';
% MaskPrefix = 'FH2_PTA_12_8_21_';
% EnergyPrefix = 'Heart2Original';
% Njm = 2048; Nim = 2048; Nkm = 1024;
% MRange = 141+[100:1123];

% Upper ventricle set
% InputPath = '../Data/STBinaryUpperVentricles/';
% OutputPath = '../Data/STBinaryUpperVentricles/';
% MaskPath = '../Images/Mask_Part2_CleanMask/AllFramesMask/';
% MaskPrefix = 'FH2_PTA_12_8_21_';
% Njm = 2048; Nim = 2048; Nkm = 512;
% MRange = [1150:1661];

% % Subregion septum
% InputPath = '../Data/STBinarySubregionSeptum/';
% OutputPath = '../Data/STBinarySubregionSeptum/';
% Njm = 150; Nim = 1122; Nkm = 794;
% MRange = [1:Nkm];
% MaskPath = '../Images/SubregionSeptum/';
% MaskPrefix = 'SubregionSeptum';
% 
% % Subregion LV freewall
% InputPath = '../Data/STBinarySubregionLVFreewall/';
% OutputPath = '../Data/STBinarySubregionLVFreewall/';
% Njm = 484; Nim = 162; Nkm = 860;
% MRange = [1:Nkm];
% MaskPath = '../Images/SubregionLVFreewall/';
% MaskPrefix = 'SubregionLVFree';

% Subregion LV freewall
% InputPath = '../Data/STBinarySubregionLVFreewallInbuiltG/';
% OutputPath = '../Data/STBinarySubregionLVFreewallInbuiltG/';
% Njm = 484; Nim = 162; Nkm = 860;
% MRange = [1:Nkm];
% MaskPath = '../Images/SubregionLVFreewall/';
% MaskPrefix = 'SubregionLVFree';

% Full heart 1
% InputPath = '../DataHeart1/STBinary/';
% OutputPath = '../DataHeart1/STBinary/';
% %HelixImageOutputPath = '../ImagesHeart1/ResizedDisplay/HelixAngleImages/Level6/';
% HelixImageOutputPath = '../ImagesHeart1/ResizedDisplay/HelixAngleImages/Level5/';
% MaskPath = '../ImagesHeart1/MaskEdited/';
% MaskPrefix = 'FH1_PTA_20_7_21_';
% EnergyPrefix = 'Heart1Original';
% Njm = 2048; Nim = 2352; Nkm = 1493;
% MRange = [260:1752];
% NewSize = [512,588,373];
% 
% Air Dry Heart 2
% InputPath = '../DataHeartDry/STBinary/';
% OutputPath = '../DataHeartDry/STBinary/';
% HelixImageOutputPath = sprintf('../ImagesHeartDry/HelixAngleImages/Level%1d/',Level);
% MaskPath = '../ImagesHeartDry/Mask/';
% MaskPrefix = 'AirDry_';
% EnergyPrefix = 'Heart2AirDry';
% Njm = 1237; Nim = 1237;
% MRange = [548:1060];
% Nkm = length(MRange);
% NewSize = [512,512,64];

% Air Dry H1C1H
% InputPath = '../Images_H1C1H/Data_H1C1H/STBinary/';
% OutputPath = '../Images_H1C1H/Data_H1C1H/STBinary/';
% HelixImageOutputPath = sprintf('../Images_H1C1H/HelixAngleImages/Level%1d/',Level);
% WriteHelixAngleImages = 'no';
% MaskPath = '../Images_H1C1H/MaskClean/';
% MaskPrefix = 'H1C1H_';
% EnergyPrefix = 'H1C1H_';
% HelixAngleFile = 'LVAxisDataForHelixAngleCalcs';
% Njm = 1589; Nim = 1771;
% MRange = [370:1689];
% Nkm = length(MRange);
% NewSize = [512,571,425];

% Air Dry H1FGR1H
InputPath = '../Images_H1FGR1H/Data_H1FGR1H/STBinary/';
OutputPath = '../Images_H1FGR1H/Data_H1FGR1H/STBinary/';
HelixImageOutputPath = sprintf('../Images_H1FGR1H/',Level);
WriteHelixAngleImages = 'yes';
MaskPath = '../Images_H1FGR1H/Mask/';
MaskPrefix = 'H1FGR1H_';
EnergyPrefix = 'H1FGR1H_';
HelixAngleFile = 'LVAxisDataForHelixAngleCalcs';
Njm = 1586; Nim = 2083;
MRange = [1:1460];
Nkm = length(MRange);
NewSize = [512,521,396];

% Air Dry H1C1H 20 Sept 2022
% InputPath = '../Images_H1C1H/Data_H1C1H/STBinary/';
% OutputPath = '../Images_H1C1H/Data_H1C1H/STBinary/';
% HelixImageOutputPath = sprintf('../Images_H1C1H/',Level);
% WriteHelixAngleImages = 'yes';
% MaskPath = '../Images_H1C1H/Mask/';
% MaskPrefix = 'H1C1H_';
% EnergyPrefix = 'H1C1H_';
% HelixAngleFile = 'LVAxisDataForHelixAngleCalcs';
% Njm = 1613; Nim = 2118;
% MRange = [1:1616];
% Nkm = length(MRange);
% NewSize = [512,530,403];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' ... Loading in data ... \n');

fid = fopen(sprintf('%sS%1d.bin',InputPath,Level),'r');
N = fread(fid,3,'uint16');
fprintf(' ... Dimension: [%d,%d,%d] ... \n',N);
d2Xs = double(fread(fid,prod(N),'double'));
dXYs = double(fread(fid,prod(N),'double'));
dXZs = double(fread(fid,prod(N),'double'));
d2Ys = double(fread(fid,prod(N),'double'));
dYZs = double(fread(fid,prod(N),'double'));
d2Zs = double(fread(fid,prod(N),'double'));
fclose(fid);

fid = fopen(sprintf('%sIJK%1d.bin',InputPath,Level),'r');
N = fread(fid,3,'uint16');
fprintf(' ... Dimension: [%d,%d,%d] ... \n',N);
I = fread(fid,prod(N),'uint16');
J = fread(fid,prod(N),'uint16');
K = fread(fid,prod(N),'uint16');
fclose(fid);

% Read in mask data
kstart = MRange(1); kend = MRange(Nkm); 
I3D = true(Njm,Nim,Nkm); 
parfor k=1:Nkm
  fnamein = sprintf('%s%s%04d.png',MaskPath,MaskPrefix,MRange(k));
  M = ~(imread(fnamein) == 0);
  I3D(:,:,k) = M;
end

% Load in precomputed axis orientations for helix angle calculations
load(sprintf('%s%s.mat',InputPath,HelixAngleFile));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Manipulate loaded data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FullMask = I3D;
I3D = permute(I3D,[2,1,3]);
I3D = reshape(I3D,Nim*Njm*Nkm,1);
GD = ((K-1)*Njm+J-1)*Nim+I;
MaskGD = find(I3D(GD));

% Only use data within masked tissue
d2Xs = d2Xs(MaskGD);
dXYs = dXYs(MaskGD);
dXZs = dXZs(MaskGD);
d2Ys = d2Ys(MaskGD);
dYZs = dYZs(MaskGD);
d2Zs = d2Zs(MaskGD);
I = I(MaskGD);
J = J(MaskGD);
K = K(MaskGD);

% Index step sizes 
DI = I(2)-I(1); DJ = J(1+N(1))-J(1); DK = K(1+N(1)*N(2))-K(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Eigenanalysis.
%
% Eigenvalues and eigenvectors are found for each tissue point (background
% points are excluded). Angles from the helix plane are found at the same
% time.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' ... Finding eigen-things ... \n');

L1 = zeros(length(d2Xs),1);
L2 = zeros(length(d2Xs),1);
L3 = zeros(length(d2Xs),1);
E1 = zeros(length(d2Xs),3);
E2 = zeros(length(d2Xs),3);
E3 = zeros(length(d2Xs),3);
HelixF = zeros(length(d2Xs),1);
HelixS = zeros(length(d2Xs),1);
HelixN = zeros(length(d2Xs),1);
TLVV = zeros(length(d2Xs),3);
rLVV = zeros(length(d2Xs),3);
parfor i=1:length(d2Xs)
  if ~mod(i,100000)
    fprintf(' entry: %d\n',i);
  end
  % local structure tensor
  S = [d2Xs(i),dXYs(i),dXZs(i);dXYs(i),d2Ys(i),dYZs(i);dXZs(i),dYZs(i),d2Zs(i)];
  [V,D] = eig(S); % evect/eval in largest to smallest
  [y,idx]=sort(diag(D));
  L1(i) = D(idx(3),idx(3));
  L2(i) = D(idx(2),idx(2));
  L3(i) = D(idx(1),idx(1)); % smallest, i.e. fibre
  E1(i,:) = V(:,idx(3))';
  E2(i,:) = V(:,idx(2))';
  E3(i,:) = V(:,idx(1))'; % smallest, i.e. fibre
  
  % Compute helix angle with respect to a plane orthogonal to the principal
  % axis of the ventricular space. 
%  PA = Coefs(:,3);
  PA = Coefs(:,1);
  DPC = [[I(i),J(i),K(i)]-COM]';
  Delta = dot(DPC,PA);
  r = DPC - Delta*PA; r = r/norm(r);
  T = cross(PA,r); T = T/norm(T);
  if dot(E3(i,:)',T) < 0 % flip the fibre direction if is is pointing toward the apex.
      E3(i,:) = -1*E3(i,:);
  end
  HelixF(i) = asind(dot(E3(i,:)',PA)/sqrt(dot(E3(i,:)',PA)^2+dot(E3(i,:)',T)^2));
  HelixS(i) = asind(dot(E2(i,:)',PA)/sqrt(dot(E2(i,:)',PA)^2+dot(E2(i,:)',T)^2));
  HelixN(i) = asind(dot(E1(i,:)',PA)/sqrt(dot(E1(i,:)',PA)^2+dot(E1(i,:)',T)^2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate information content using Shannon entropy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' ... Calculate information content ... \n');

% Calculate a fractional anisotropy - for Shannon entropy calculation
Trace = (L1+L2+L3)/3;
IG = (L1-L3 <= 1e-6) & (L1+L2+L3 <= 1e-6);
Denom = sqrt(L1.^2+L2.^2+L3.^3+1e-6*cast(IG,'double'));
FA1 = sqrt(3/2)*(sqrt((L1-Trace).^2+(L2-Trace).^2+(L3-Trace).^2))./Denom;

% Structure tensor data
STD = [d2Xs,dXYs,dXZs,d2Ys,dYZs,d2Zs];

% Compute information content - probability histograms first and then
% compute Shannon entropy. Shannon entropy is the amount of information in
% a variable.
% Probability distributions
InfoLabels = {'Eigenvalues';'Eigenvectors';'Fibre Helix Angles';'Fractional Anisotropy';'Structure Tensor'};
HEval = histogram([L1,L2,L3],1000,'Normalization','Probability'); HEval = HEval.Values(find(abs(HEval.Values) > 0));
HEvec = histogram([E1,E2,E3],1000,'Normalization','Probability'); HEvec = HEvec.Values(find(abs(HEvec.Values) > 0));
HHF = histogram(HelixF,1000,'Normalization','Probability'); HHF = HHF.Values(find(abs(HHF.Values) > 0));
HFA = histogram(FA1,1000,'Normalization','Probability'); HFA = HFA.Values(find(abs(HFA.Values) > 0));
HSTD = histogram(STD,1000,'Normalization','Probability'); HSTD = HSTD.Values(find(abs(HSTD.Values) > 0));
% Compute Shannon entropy
InfoContent = [-HEval*log2(HEval');-HEvec*log2(HEvec');-HHF*log2(HHF');-HFA*log2(HFA');-HSTD*log2(HSTD')];
% Save Shannon entropy
save(sprintf('%s%s_ShannonLogEnergyEntropyData_L%d.mat',OutputPath,EnergyPrefix,Level),'InfoLabels','InfoContent');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Assess eigenvalues and intensity gradients
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' ... Assessing eigenvalues and gradients ... \n');

% compute measures of strength of eigen values
CO = zeros(length(d2Xs),1);
COFS = zeros(length(d2Xs),1);
COFN = zeros(length(d2Xs),1);
COSN = zeros(length(d2Xs),1);
IG = zeros(length(d2Xs),1);
CL = zeros(length(d2Xs),1);
GN = zeros(length(d2Xs),1); % norm of the gradient

IG = (L1-L3 <= 1e-6) & (L1+L2+L3 <= 1e-6);

COFS = (L2-L3)./max(L2+L3,1e-6);
COFN = (L1-L3)./max(L1+L3,1e-6);
COSN = (L1-L2)./max(L1+L2,1e-6);

CO = (COFS+COFN+COSN)/3;
CON = (COFN+COSN)/2;

% Calculate a coherency index
CO = (L1-L3)./(L1+L2+L3+1e6*cast(IG,'double'));
CL = ((L1>1e-6 & L2>1e-6 & L3>1e-6)*3) + ((L1>1e-6 & L2>1e-6 & L3<=1e-6)*2) + ((L1>1e-6 & L2<=1e-6 & L3<=1e-6)*1);

% Calculate a fractional anisotropy
Trace = (L1+L2+L3)/3;
Denom = sqrt(L1.^2+L2.^2+L3.^3+1e-6*cast(IG,'double'));
FA1 = sqrt(3/2)*(sqrt((L1-Trace).^2+(L2-Trace).^2+(L3-Trace).^2))./Denom;
FA2 = sqrt(1/2)*(sqrt((L1-L2).^2+(L2-L3).^2+(L3-L1).^2))./Denom;

% Norm of the intensity gradient
GN = sqrt(d2Xs+d2Ys+d2Zs);
IFG = (sum([sqrt(abs(d2Xs)),sqrt(abs(d2Ys)),sqrt(abs(d2Zs))].*E3,2));


ttotalprocess1 = etime(clock,ttotalprocess0); fprintf(' *** total processing time: %0.2f sec\n',ttotalprocess1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Write text file in exdata format for display in cmgui
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' ... Writing display exdata file ... \n');
Idx = 1:length(I);
NLayer = [N(1),N(2),N(3)];

% Organise data
Is = I(Idx);
Js = J(Idx);
Ks = K(Idx);
Xs = Is; 
Ys = Js;
Zs = Ks;
E1s = E1(Idx,:);
E2s = E2(Idx,:);
E3s = E3(Idx,:);
COs = CO(Idx);
COFSs = COFS(Idx);
COFNs = COFN(Idx);
COSNs = COSN(Idx);
L1s = L1(Idx);
L2s = L2(Idx);
L3s = L3(Idx);
FA1s = FA1(Idx);
FA2s = FA2(Idx);
Angles = HelixF(Idx);
GNorm = GN(Idx);
IntensityFGradient = IFG(Idx);
HelixFs = HelixF(Idx);
HelixSs = HelixS(Idx);
HelixNs = HelixN(Idx);

% Write file
WriteExdataFile(NLayer,1,Xs,Ys,Zs,E1s,E2s,E3s,COs,COFSs,COFNs,COSNs,L1s,L2s,L3s,FA1s,FA2s,Angles,GNorm,IntensityFGradient,HelixFs,HelixSs,HelixNs,sprintf('%sVisData_L%d.exdata',OutputPath,Level),'OrientationF',300000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determine helix angles in an image format for display purposes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('... Interpolate helix angles onto a resized image ... \n');

% Resize mask images
RFullMask = imresize3(FullMask,NewSize);
% Sample points for helix angle image
[RFMJ,RFMI,RFMK]=ndgrid(round(1:(Njm/NewSize(1)):Njm),round(1:(Nim/NewSize(2)):Nim),round(1:(Nkm/NewSize(3)):Nkm));
% Construct interpolating function for fibre helix angle
FSI = scatteredInterpolant(I,J,K,HelixF,'natural','linear');

% Write interpolated image stack
if strcmp(WriteHelixAngleImages,'yes')
    fprintf('... Writing interpolated images ...\n')
    parfor k=1:NewSize(3)
      if ~mod(k,20) fprintf(' slice: %d\n',k); end
      Slice = FSI(RFMI(:,:,k),RFMJ(:,:,k),RFMK(:,:,k));
      ImSlice = uint8(255*(Slice+90)/180).*uint8(RFullMask(:,:,k));
      fname = sprintf('%sHA%03d.png',HelixImageOutputPath,k);
      imwrite(ImSlice,fname);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Export layer of points around the LV COM.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find layer number nearest to COM_k, i.e. COM(3)
UK = unique(K);
DCOM = abs(UK-repmat(round(COM(3)),[length(UK),1]));
[MV,VMIdx] = min(DCOM);
Idx = find(K==UK(VMIdx));
exfname = sprintf('%sVisLVCOMLayerData_L%d',OutputPath,Level);
DataSLabels = {'FA1','HelixF'};
DataVLabels = {'Fiber'};
DataS = zeros(length(Idx),2);
DataS(:,1) = FA1(Idx);
DataS(:,2) = HelixF(Idx);
DataV = cell(1);
DataV{1} = E3(Idx,:);
GName = sprintf('HelixVector');
WriteGeneralExdataFile(I(Idx),J(Idx),K(Idx),[1:length(Idx)]',DataS,DataV,exfname,GName,DataSLabels,DataVLabels);




