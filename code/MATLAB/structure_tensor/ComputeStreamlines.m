%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ComputeStreamlines.m
%
% This script loads and processes structure tensor components at a 
% specified resolution to determine streamlines as a visible expression of 
% fibre direction.
%
% Updated by: Mark Trew, June 2022.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

ttotalprocess0 = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set parameters and paths
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Level = 5; % frequency resolution of ST/Hessian data to use

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Data and path locations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
base_dir = join([getenv("HOME"), "Documents/phd"], '/');
InputPath = base_dir + '/AWA015_PTA_1_Rec_Trans/downsampled/ST/binary/';
OutputPath = base_dir + '/AWA015_PTA_1_Rec_Trans/downsampled/ST/binary/';
MaskPath = base_dir + '/AWA015_PTA_1_Rec_Trans/downsampled/ST/mask/';
MaskPrefix = 'AWA015_PTA_1_';
kindex = 80:534; 
Nj = 548; Ni = 1124;
Nkm = length(kindex);
%%
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% InputPath = base_dir + '/AWA015_PTA_2_Ova_Rec_Trans/downsampled/ST/binary/';
% OutputPath = base_dir + '/AWA015_PTA_2_Ova_Rec_Trans/downsampled/ST/binary/';
% MaskPath = base_dir + '/AWA015_PTA_2_Ova_Rec_Trans/downsampled/ST/mask/';
% MaskPrefix = 'AWA015_PTA_2_Ova_Trans_';
% kindex = 1:325; 
% Nj = 518; Ni = 753;
% Nkm = length(kindex);
%% Downsampled version
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% InputPath = base_dir + '/AWA014_PTA_2_Rec_Trans/downsampled/ST/binary/';
% OutputPath = base_dir + '/AWA014_PTA_2_Rec_Trans/downsampled/ST/binary/';
% MaskPath = base_dir + '/AWA014_PTA_2_Rec_Trans/downsampled/ST/mask/';
% MaskPrefix = 'AWA014_PTA_2_';
% Nj = 469; Ni = 512;
% kindex = 1:138;
% Nkm = length(kindex);

%% Full resolution segment version 
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% InputPath = base_dir + '/AWA014_PTA_2_Rec_Trans/ST/binary/';
% OutputPath = base_dir + '/AWA014_PTA_2_Rec_Trans/ST/binary/';
% MaskPath = base_dir + '/AWA014_PTA_2_Rec_Trans/ST/mask/';
% MaskPrefix = 'AWA014_PTA_2_';
% kindex = 830:1030; 
% Nj = 1401; Ni = 1001;
% Nkm = length(kindex);

%% Reoriented downsampled version 
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% InputPath = base_dir + '/AWA014_PTA_2_Rec_Trans/downsampled/ST/binary/';
% OutputPath = base_dir + '/AWA014_PTA_2_Rec_Trans/downsampled/ST/binary/';
% MaskPath = base_dir + '/AWA014_PTA_2_Rec_Trans/downsampled/ST/mask/';
% MaskPrefix = 'AWA014_PTA_2_';
% Nj = 367; Ni = 473;
% kindex = 1:372;
% Nkm = length(kindex);

%% Reoriented full resolution segment version 
% base_dir = join([getenv("HOME"), "Documents/phd"], '/');
% InputPath = base_dir + '/AWA014_PTA_2_Rec_Trans/ST/binary/';
% OutputPath = base_dir + '/AWA014_PTA_2_Rec_Trans/ST/binary/';
% MaskPath = base_dir + '/AWA014_PTA_2_Rec_Trans/ST/mask/';
% MaskPrefix = 'AWA014_PTA_2_';
% Nj = 420; Ni = 388;
% kindex = 1:383;
% Nkm = length(kindex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' ... loading in data ... \n');

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

% Read in mask data and 
kstart = kindex(1); kend = kindex(Nkm); 
I3D = true(Nj,Ni,Nkm); 
parfor k=1:Nkm
  fnamein = sprintf('%s%s%03d.png',MaskPath,MaskPrefix,kindex(k));
  M = ~(imread(fnamein) == 0);
  I3D(:,:,k) = M;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Manipulate loaded data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FullMask = I3D;
I3D = permute(I3D,[2,1,3]);
I3D = reshape(I3D,Ni*Nj*Nkm,1);
GD = ((K-1)*Nj+J-1)*Ni+I;
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
DJ = J(1+N(1))-J(1); DK = 225; %DK = K(1+N(1)*N(2))-K(1);
DI = 32; % fix to this value regardless of data level.

fprintf('DI: %d, DJ: %d, DK: %d\n',DI,DJ,DK);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set up sample grid (seed points) for streamlines. The density of the seed
% points is arbitrary. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('... Set up sample grid: ');

[SI,SJ,SK] = ndgrid(1:round(2*DI/0.5):Ni,1:round(2*DI/0.5):Nj,1:round(2*DI/0.5):Nkm);
%[SI,SJ,SK] = ndgrid(1:round(2*DI/1):Nim,1:round(2*DI/1):Njm,1:round(2*DI/1):Nkm);
%[SI,SJ,SK] = ndgrid(1:round(2*DI/2):Nim,1:round(2*DI/2):Njm,1:round(2*DI/2):Nkm);
fprintf('%dX%dX%d\n',size(SI));

NS = size(SI);
LIS = sub2ind([Ni,Nj,Nkm],reshape(SI,[prod(NS),1]),reshape(SJ,[prod(NS),1]),reshape(SK,[prod(NS),1]));
IdxS = find(I3D(LIS));

% Set up interpolant
Fd2Xs = scatteredInterpolant(I,J,K,d2Xs);
FdXYs = scatteredInterpolant(I,J,K,dXYs);
FdXZs = scatteredInterpolant(I,J,K,dXZs);
Fd2Ys = scatteredInterpolant(I,J,K,d2Ys);
FdYZs = scatteredInterpolant(I,J,K,dYZs);
Fd2Zs = scatteredInterpolant(I,J,K,d2Zs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate streamline paths
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine paths
Paths = cell(length(IdxS),1);
DS = 5;

%FiberIndex = 2; % sheet - otherwise use 1 for smallest eigenvalue/fiber
%FiberIndex = 3; % normal - otherwise use 1 for smallest eigenvalue/fiber
FiberIndex = 1; % fiber - otherwise use 1 for smallest eigenvalue/fiber
MaxTrackLength = 10000; % Fiber tracks
%MaxTrackLength = 500; % sheet tracks
parfor i=1:length(IdxS)
%for i=1:1
    if ~mod(i,10) fprintf('Path: %d\n',i); end
    Paths{i} = FiberTrack([SI(IdxS(i)),SJ(IdxS(i)),SK(IdxS(i))],DS,I,J,K,Fd2Xs,FdXYs,FdXZs,Fd2Ys,FdYZs,Fd2Zs,I3D,[Ni,Nj,Nkm],FiberIndex,MaxTrackLength);
%                    L3 = Paths{i}{1}(n,1); L2 = Paths{p}{d}(n,2); L1 = Paths{p}{d}(n,3);
%                   Trace = (L1+L2+L3)/3;
%                   Denom = sqrt(L1.^2+L2.^2+L3.^3+1e-8);
%                   FA1 = sqrt(3/2)*(sqrt((L1-Trace).^2+(L2-Trace).^2+(L3-Trace).^2))./Denom;
end

% % Try and orient paths  
% EndPoints = zeros(length(Paths),4,2);
% for i=1:length(Paths)
%   LP = [flipud(Paths{i}{2});Paths{i}{1}];
%   EndPoints(i,1:4,1) = [LP(1,1),LP(1,2),LP(1,3),0];
%   EndPoints(i,1:4,2) = [LP(end,1),LP(end,2),LP(end,3),1];
% end

%exfname = sprintf('%s/Streamlines',OutputPath);
if FiberIndex == 2
  exfname = sprintf('%s/StreamlinesSheet',OutputPath);
elseif FiberIndex == 3
  exfname = sprintf('%s/StreamlinesNormal',OutputPath);
else
  exfname = sprintf('%s/Streamlines_L%1d_FB',OutputPath,Level);
end
groupname = 'Streamlines';
ExportStreamlines(Paths,[],exfname,groupname,1,1);

