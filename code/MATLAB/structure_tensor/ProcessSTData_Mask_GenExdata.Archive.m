% This script loads in and processes structure tensor components 
% at a specifield frequency level.

clear all;

ttotalprocess0 = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set parameters and paths
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Level = 5; % frequency resolution of ST/Hessian data to use

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InputPath = '../Data/STBinary/';
OutputPath = '../Data/STBinary/';
%MaskPath = '../RawImages/BinaryMask/';
MaskPath = '../Images/Mask/';
MaskPrefix = 'FH2_PTA_12_8_21_';
Njm = 2048; Nim = 2048; Nkm = 1024;
MRange = 141+[100:1123];

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
kstart = MRange(1); kend = MRange(Nkm); 
I3D = true(Njm,Nim,Nkm); 
for k=kstart:kend,
  fnamein = sprintf('%s%s%04d.png',MaskPath,MaskPrefix,k);
  M = (imread(fnamein)==1);
  %%M = imerode(M,ones(13,13));
  I3D(:,:,k-kstart+1) = M;
end;
FullMask = I3D;
I3D = permute(I3D,[2,1,3]);
Distance = bwdist(~I3D);
I3D = reshape(I3D,Nim*Njm*Nkm,1);
GD = ((K-1)*Njm+J-1)*Nim+I;
MaskGD = find(I3D(GD));

%I3D = true(Njm,Nim,Nkm); 
%I3D = reshape(I3D,Nim*Njm*Nkm,1);
%GD = ((K-1)*Njm+J-1)*Nim+I;
%MaskGD = find(I3D(GD));

% Only use data within masked tissue
d2Xs = d2Xs(MaskGD);
dXYs = dXYs(MaskGD);
dXZs = dXZs(MaskGD);
d2Ys = d2Ys(MaskGD);
dYZs = dYZs(MaskGD);
d2Zs = d2Zs(MaskGD);
dXs = dXs(MaskGD);
dYs = dYs(MaskGD);
dZs = dZs(MaskGD);
%h2Xs = h2Xs(MaskGD);
%hXYs = hXYs(MaskGD);
%hXZs = hXZs(MaskGD);
%h2Ys = h2Ys(MaskGD);
%hYZs = hYZs(MaskGD);
%h2Zs = h2Zs(MaskGD);

% build adjacency links so that 
%IOrig = I; JOrig = J; KOrig = K;
DI = I(2)-I(1); DJ = J(1+N(1))-J(1); DK = K(1+N(1)*N(2))-K(1);

I = I(MaskGD);
J = J(MaskGD);
K = K(MaskGD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Sample grid for streamlines
% %[SI,SJ,SK] = ndgrid(1:round(DI/3):Nim,1:round(DI/3):Njm,1:round(DI/3):Nkm);
% %[SI,SJ,SK] = ndgrid(1:round(DI/2):Nim,1:round(DI/2):Njm,1:round(DI/2):Nkm);
% %[SI,SJ,SK] = ndgrid(1:round(DI/1):Nim,1:round(DI/1):Njm,1:round(DI/1):Nkm);
% [SI,SJ,SK] = ndgrid(1:2*round(DI/1):Nim,1:2*round(DI/1):Njm,1:2*round(DI/1):Nkm);
% % Find linear indices
% NS = size(SI);
% LIS = sub2ind([Nim,Njm,Nkm],reshape(SI,[prod(NS),1]),reshape(SJ,[prod(NS),1]),reshape(SK,[prod(NS),1]));
% IdxS = find(I3D(LIS));
% 
% % Set up interpolant
% Fd2Xs = scatteredInterpolant(I,J,K,d2Xs);
% FdXYs = scatteredInterpolant(I,J,K,dXYs);
% FdXZs = scatteredInterpolant(I,J,K,dXZs);
% Fd2Ys = scatteredInterpolant(I,J,K,d2Ys);
% FdYZs = scatteredInterpolant(I,J,K,dYZs);
% Fd2Zs = scatteredInterpolant(I,J,K,d2Zs);
% 
% % Determine paths
% Paths = cell(length(IdxS),1);
% DS = 5;
% parfor i=1:length(IdxS)
% %for i=1:1
%   if ~mod(i,100) fprintf('Path: %d\n',i); end
%   Paths{i} = FiberTrack([SI(IdxS(i)),SJ(IdxS(i)),SK(IdxS(i))],DS,I,J,K,Fd2Xs,FdXYs,FdXZs,Fd2Ys,FdYZs,Fd2Zs,I3D,[Nim,Njm,Nkm]);
% end
% 
% % Try and orient paths  
% EndPoints = zeros(length(Paths),4,2);
% for i=1:length(Paths)
%   LP = [flipud(Paths{i}{2});Paths{i}{1}];
%   EndPoints(i,1:4,1) = [LP(1,1),LP(1,2),LP(1,3),0];
%   EndPoints(i,1:4,2) = [LP(end,1),LP(end,2),LP(end,3),1];
% end
% 
% % Loop over data to cluster start points
% TestPaths = [1:size(EndPoints,1)]';
% CurrentPoints = [1:size(EndPoints,1)]';
% UnflippedStartPoints = [];
% FlippedStartPoints = [];
% while ~isempty(TestPaths)
%     fprintf('Length test paths: %d\n',length(TestPaths));
%     ClusterEndPoints = [EndPoints(TestPaths,1:3,1);EndPoints(TestPaths,1:3,2)];
%     CIdx = kmeans(ClusterEndPoints,2);
%     [dummy,MinSizeCluster] = min([length(find(CIdx==1)),length(find(CIdx==2))]);
%     IdxTCluster = find(CIdx==MinSizeCluster);
%     WorkingUnflipped = IdxTCluster(find(IdxTCluster <= size(TestPaths,1)));
%     WorkingFlipped = IdxTCluster(find(IdxTCluster > size(TestPaths,1)))-size(TestPaths,1);
%     WorkingFlipped = setdiff(WorkingFlipped,WorkingUnflipped);
%     TestPaths = setdiff([1:size(TestPaths,1)],[WorkingUnflipped;WorkingFlipped])';
%     UnflippedStartPoints = [UnflippedStartPoints;CurrentPoints(WorkingUnflipped)];
%     FlippedStartPoints = [FlippedStartPoints;CurrentPoints(WorkingFlipped)];
%     CurrentPoints = CurrentPoints(TestPaths);
% end
% 
% % Update endpoints
% UEndPoints = EndPoints;
% for i=1:length(FlippedStartPoints)
%     Temp = UEndPoints(FlippedStartPoints(i),:,1);
%     UEndPoints(FlippedStartPoints(i),:,1) = [UEndPoints(FlippedStartPoints(i),1:3,2),0];
%     UEndPoints(FlippedStartPoints(i),:,2) = [Temp(1,1:3),1];
% end
% FullEndPoints = [UEndPoints(:,:,1);UEndPoints(:,:,2)];
% %AdjS = cell(size(FullEndPoints,1),1);
% AdjS = cell(size(UEndPoints,1),1);
% for i=1:length(AdjS)
%   D = sum([(FullEndPoints(:,1)-UEndPoints(i,1,2)),(FullEndPoints(:,2)-UEndPoints(i,2,2)),(FullEndPoints(:,3)-UEndPoints(i,3,2))].^2,2);
%   INearest = find(D <= (3*DS)^2);
%   AdjS{i} = INearest;
% end
% % find where there are major differences, i.e. end points in mainly start
% % point regions and vice versa
% CompareMedian = zeros(size(AdjS,1),1);
% for i=1:size(AdjS,1)
%   CompareMedian(i) = median(FullEndPoints(AdjS{i},4));
% end
% NewFlippedStartPoints = union(FlippedStartPoints,find(CompareMedian==0));
% NewUnflippedStartPoints = setdiff(UnflippedStartPoints,find(CompareMedian==0));

% Draw streamlines
% DrawStreamlines(Paths,NewUnflippedStartPoints,NewFlippedStartPoints,
% DrawStreamlines(Paths,NewUnflippedStartPoints,NewFlippedStartPoints,Distance,2,[20 100 20 80 20 120]);
% DrawStreamlines(Paths,NewUnflippedStartPoints,NewFlippedStartPoints,Distance,2,[]);






% Notes: streamlines are the oriented fiber field. Next stage is to project
% the streamlines back to the grid (or something similar). Need to account for
% the flippped fiber vectors.
% A flipped fiber vector should have a corresponding flipped sheet vector.
% Possibly as well as returning the path, the FiberTrack should also return 
% fiber, sheet and normal vectors with sheet flipped when fiber is flipped.
% Can then also apply 
% Fibers, sheets and normal vectors can then also be taken back to regular
% (I,J,K) grid (with mask) using a nearest interpolant etc.
%
% Problem is that while consistent fiber orientations along a streamline 
% can be found, these may not be consistent from streamline to streamline,
% so whole thing is back to the same issue.
%
% Clustering approach help a little, but not entirely satisfactory. Perhaps 
% another option 

%return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Eigenanalysis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' ... finding ethings ... \n');
L1 = zeros(length(d2Xs),1);
L2 = zeros(length(d2Xs),1);
L3 = zeros(length(d2Xs),1);
E1 = zeros(length(d2Xs),3);
E2 = zeros(length(d2Xs),3);
E3 = zeros(length(d2Xs),3);
L1H = zeros(length(h2Xs),1);
L2H = zeros(length(h2Xs),1);
L3H = zeros(length(h2Xs),1);
E1H = zeros(length(h2Xs),3);
E2H = zeros(length(h2Xs),3);
E3H = zeros(length(h2Xs),3);
parfor i=1:length(d2Xs),
  if ~mod(i,100000)
    fprintf(' entry: %d\n',i);
  end;
  % local structure tensor
  S = [d2Xs(i),dXYs(i),dXZs(i);dXYs(i),d2Ys(i),dYZs(i);dXZs(i),dYZs(i),d2Zs(i)];
  [V,D] = eig(S); % evect/eval in largest to smallest
  [y,idx]=sort(diag(D));
  L1(i) = D(idx(3),idx(3));
  L2(i) = D(idx(2),idx(2));
  L3(i) = D(idx(1),idx(1));
  E1(i,:) = V(:,idx(3))';
  E2(i,:) = V(:,idx(2))';
  E3(i,:) = V(:,idx(1))';

%  H = [h2Xs(i),hXYs(i),hXZs(i);hXYs(i),h2Ys(i),hYZs(i);hXZs(i),hYZs(i),h2Zs(i)];
%  [VH,DH] = eig(H); % evect/eval in largest to smallest
%  [yH,idxH]=sort(diag(DH));
%  L1H(i) = DH(idxH(3),idxH(3));
%  L2H(i) = DH(idxH(2),idxH(2));
%  L3H(i) = DH(idxH(1),idxH(1));
%  E1H(i,:) = VH(:,idxH(3))';
%  E2H(i,:) = VH(:,idxH(2))';
%  E3H(i,:) = VH(:,idxH(1))';
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Assess eigenvalues
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coherency and categorisation
fprintf(' ... categorizing ... \n');

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

% Path track
%[Path]=FiberTrack([410,275,325],1,100,I,J,K,d2Xs,dXYs,dXZs,d2Ys,dYZs,d2Zs);

%
TestAngle = (180/pi)*asin(E3(:,3)./sqrt(sum(E3(:,[1,3]).^2,2)));
Q1 = double(E3(:,3) >=0 & E3(:,1) >0);
Q2 = double(E3(:,3) >=0 & E3(:,1) <= 0);
Q3 = double(E3(:,3) <0  & E3(:,1) < 0);
Q4 = double(E3(:,3) <0  & E3(:,1) >=0);
FiberAngles = Q1.*TestAngle + Q2.*(-TestAngle) + Q3.*(-TestAngle) + Q4.*(TestAngle);
FiberAngles = -1*FiberAngles; % Convention

% Flip Fiber vectors as necessary:
%E3 = repmat(Q1,[1,3]).*E3 + repmat(Q2,[1,3]).*(-1*E3) +repmat(Q3,[1,3]).*(-1*E3) + repmat(Q4,[1,3]).*(E3);

% absolute intensity gradient magnitude in the fiber direction
IFG = (sum([dXs,dYs,dZs].*E3,2));
ISG = (sum([dXs,dYs,dZs].*E2,2));
ING = (sum([dXs,dYs,dZs].*E1,2));

ttotalprocess1 = etime(clock,ttotalprocess0); fprintf(' *** total processing time: %0.2f sec\n',ttotalprocess1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Write files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filename = sprintf('%sFiber_L%d.dat',OutputPath,Level);
%fid = fopen(filename,'w');
%fprintf(fid,'%f %f %f\n',E3');
%fclose(fid);

%filename = sprintf('%sSheet_L%d.dat',OutputPath,Level);
%fid = fopen(filename,'w');
%fprintf(fid,'%f %f %f\n',E2');
%fclose(fid);

%filename = sprintf('%sNormal_L%d.dat',OutputPath,Level);
%fid = fopen(filename,'w');
%fprintf(fid,'%f %f %f\n',E1');
%fclose(fid);

%filename = sprintf('%sEvalues_L%d.dat',OutputPath,Level);
%fid = fopen(filename,'w');
%fprintf(fid,'%f %f %f\n',[L3,L2,L1]');
%fclose(fid);

% Save matlab data
save(sprintf('%sData_L%d.mat',OutputPath,Level),'I','J','K','E1','E2','E3','L1','L2','L3','IG','COFS','COFN','COSN','CO','CL','Trace','FA1','GN','FiberAngles','IFG','ISG','ING','-v7.3');


%[mv,Idx] = min(abs(K-median(unique(K))));
%Idx = find(K==K(Idx));
%NLayer = [N(1),N(2),1];



Idx = 1:length(I);
NLayer = [N(1),N(2),N(3)];

Is = I(Idx);
Js = J(Idx);
Ks = K(Idx);
Xs = Is; % modify to represent mapping from indices to xyz
Ys = Js;
Zs = Ks;
E1s = E1(Idx,:);
E2s = E2(Idx,:);
E3s = E3(Idx,:);
% E1Os = E1_Orig(Idx,:);
% E2Os = E2_Orig(Idx,:);
% E3Os = E3_Orig(Idx,:);
%E1Hs = E1H(Idx,:);
%E2Hs = E2H(Idx,:);
%E3Hs = E3H(Idx,:);
COs = CO(Idx);
COFSs = COFS(Idx);
COFNs = COFN(Idx);
COSNs = COSN(Idx);
L1s = L1(Idx);
L2s = L2(Idx);
L3s = L3(Idx);
%L1Hs = L1H(Idx);
%L2Hs = L2H(Idx);
%L3Hs = L3H(Idx);
FA1s = FA1(Idx);
FA2s = FA2(Idx);
Angles = FiberAngles(Idx);
GNorm = GN(Idx);
%GNorm = WB(Idx);
%GNorm = OldLocalObj(Idx);
IntensityFGradient = IFG(Idx);

% Turn off writing exdata file
%WriteExdataFile(NLayer,1,Xs,Ys,Zs,E1Os,E2Os,E3Os,COs,COFSs,COFNs,COSNs,L1s,L2s,L3s,FA1s,FA2s,Angles,GNorm,IntensityFGradient,sprintf('%sVisData_Orig_L%d.exdata',OutputPath,Level),'OrientationS',200000);
WriteExdataFile(NLayer,1,Xs,Ys,Zs,E1s,E2s,E3s,COs,COFSs,COFNs,COSNs,L1s,L2s,L3s,FA1s,FA2s,Angles,GNorm,IntensityFGradient,sprintf('%sVisData_L%d.exdata',OutputPath,Level),'OrientationF',300000);
%GNorm = NewLocalObj(Idx);
%WriteExdataFile(NLayer,1,Xs,Ys,Zs,E1s,NewSheets,NewFibers,COs,COFSs,COFNs,COSNs,L1s,L2s,L3s,FA1s,FA2s,Angles,GNorm,IntensityFGradient,sprintf('%sVisData_Flip_L%d.exdata',OutputPath,Level),'FlipF',400000);
%WriteExdataFile(NLayer,1,Xs,Ys,Zs,E1Hs,E2Hs,E3Hs,COs,COFSs,COFNs,COSNs,L1Hs,L2Hs,L3Hs,FA1s,FA2s,Angles,GNorm,IntensityFGradient,sprintf('%sVisData_H_L%d.exdata',OutputPath,Level),'OrientationH',300000);
%WriteExdataFile(NLayer,1,Xs,Ys,Zs,[dXs,dYs,dZs],E2s,E3s,COs,COFSs,COFNs,COSNs,L1Hs,L2Hs,L3Hs,FA1s,FA2s,Angles,GNorm,IntensityFGradient,sprintf('%sVisData_G_L%d.exdata',OutputPath,Level),'OrientationG',300000);

%return;

% Find discretisation deltas in each direction
MF = false(Njm,Nim,Nkm);
UIs = unique(Is);
UJs = unique(Js);
UKs = unique(Ks);
DI = UIs(2)-UIs(1);
DJ = UJs(2)-UJs(1);
DK = UKs(2)-UKs(1);
for i=1:length(Js)
    [P,C] = FSNEllipse(L3(i),L2(i),L1(i),E3(i,:)',E2(i,:)',E1(i,:)',DJ,DI,DK);
    [JJ,II,KK] = ndgrid((Js(i)-(C(1)-1):Js(i)+(DJ-C(1))),(Is(i)-(C(2)-1):Is(i)+(DI-C(2))),(Ks(i)-(C(3)-1):Ks(i)+(DK-C(3))));
    PIdx = find(JJ >= 1 & JJ <= Njm & II >= 1 & II <= Nim & KK >= 1 & KK <= Nkm);
    GIdx = sub2ind([Njm,Nim,Nkm],JJ(PIdx),II(PIdx),KK(PIdx));
    MF(GIdx) = P(PIdx);
end


% Loop over sample points and for each find the local FSNPlane and place it
% in a mask
MF = false(Njm,Nim,Nkm);
total = 0;
for i=1:length(Js)
    [P,DF,DS] = FSNPlane(L3(i),L2(i),L1(i),E3(i,:)',E2(i,:)',E1(i,:)',max([Js(2)-Js(1),Is(2)-Is(1),Ks(2)-Ks(1)]));
%    P = imdilate(P,ones(3,3,3));
    P = imdilate(P,ones(9,9,9));
    NP = size(P); if length(NP) < 3 NP = [NP,0]; end 
    total = total+sum(P(:));
    DNP = floor((NP-1)/2); % NP is odd, so this will be even
    dummy = false(length(max(1,Js(i)-DNP(1)):min(Njm,Js(i)+DNP(1))),length(max(1,Is(i)-DNP(2)):min(Nim,Is(i)+DNP(2))),length(max(1,Ks(i)-DNP(3)):min(Nkm,Ks(i)+DNP(3))));
    MF(max(1,Js(i)-DNP(1)):min(Njm,Js(i)+DNP(1)),max(1,Is(i)-DNP(2)):min(Nim,Is(i)+DNP(2)),max(1,Ks(i)-DNP(3)):min(Nkm,Ks(i)+DNP(3))) = P(1:size(dummy,1),1:size(dummy,2),1:size(dummy,3));
end

IdxMF = find(MF);
[II,JJ,KK]=ind2sub(size(MF),IdxMF);
WriteGeneralExdataFile(II,JJ,KK,1:length(II),[],sprintf('%sFSNSheets_L%d',OutputPath,Level),'Sheets','')

% This approach to identifying volumes works well
%MF2 = xor(MF,imerode(~FullMask,ones(5,5,5)));
%MF2 = xor(MF,~FullMask);
%DMF = bwdist(MF2);
%DMF(~FullMask) = 0;
DMF = bwdist(MF);
DMF8 = uint8(256*(DMF-min(DMF(:)))./(max(DMF(:))-min(DMF(:))));
WSEdge = watershed(DMF8,26);
WSEdge(~FullMask) = 0;
WSEdge(~(imdilate(FullMask,ones(3,3,3)))) = max(WSEdge(:))+1;
SWSEdge = WSEdge == 0;
%SWSEdge = SWSEdge.*FullMask;

% For each edge point, identify whether it is face, edge or vertex or
% external
% External has only one adjacent non zero value
% Face as two adjacent non-zero values
% Edge has more than two, as does vertext (not sure of vertex test at the
% moment)
% First pass, just identify points with one adjacent value and two adjacent
% values.
IdxEdge = find(SWSEdge);
NewLabel = uint8(false(size(SWSEdge)));
LabelFlags = uint8(false(length(IdxEdge),1));
parfor i=1:length(IdxEdge)
    [jj,ii,kk]=ind2sub(size(SWSEdge),IdxEdge(i));
    T = WSEdge(max(1,jj-1):min(Njm,jj+1),max(1,ii-1):min(Nim,ii+1),max(1,kk-1):min(Nkm,kk+1));
    if length(unique(T)) == 1 % nothing
    elseif length(unique(T)) == 2 % external face
        LabelFlags(i) = 85;
    elseif length(unique(T)) == 3; % internal face
        LabelFlags(i) = 170;
    else
        LabelFlags(i) = 255;
    end
end
NewLabel(IdxEdge) = LabelFlags;
cmap = colormap('jet');

parfor k=1:Nkm
  imwrite(SWSEdge(:,:,k),sprintf('%s%03d.png','../RawImages/TessExampleSegment/SWS_',k),'BitDepth',1);%
%  imwrite(label2rgb(NewLabel(:,:,k),cmap),sprintf('%s%03d.png','../RawImages/TessExampleSegmentColour/SWS_',k));
end

return;

% % Fiber mask. Point at each fiber location and on pixels in +/- fiber direction. Number of pixels are 
% PMDist = round(ceil(L1./L3)/2);
% MF = false(Njm,Nim,Nkm);
% %MFIdx = sub2ind(size(MF),Js,Is,Ks));
% for i=1:length(Js)
%   for j=-PMDist(i):PMDist(i)
%     OnPoints = [max(1,min(Njm,round(Js(i)+E3(i,1)*j))),max(1,min(Nim,round(Is(i)+E3(i,2)*j))),max(1,min(Nkm,round(Ks(i)+E3(i,3)*j)))];
%     MF(OnPoints(1),OnPoints(2),OnPoints(3)) = 1;
%   end
% end
% 
% % Set up a flat plane in a small volume and then rotate in 3D
% DF = ceil(L1./L3); DS = ceil(L1./L2);



% % Load in foreground mask
% ForegroundMaskTemplate = '../RawImages/BinaryForeground/h03_';
% Ni = 756; Nj=1196; Nk=500;
% IFM = false(Ni,Nj,Nk);
% parfor k=1:Nk
%     IFM(:,:,k) = imread(sprintf('%s%03d.png',ForegroundMaskTemplate,k));
% end

DMF = bwdist(MF);
DMF = -DMF;
DMF(~FullMask) = Inf;
TessEdge = watershed(DMF);
%TessEdge(~FullMask) = 0;

% Write colormap
tcmap = colormap('parula');
DMF8 = uint8(256*(DMF-min(DMF(:)))./(max(DMF(:))-min(DMF(:))));
parfor k=1:Nkm
    Irgb = ind2rgb(DMF8(:,:,k),tcmap);
    imwrite(Irgb,sprintf('%s%03d.png','../RawImages/TessExample/',k)); 
end

