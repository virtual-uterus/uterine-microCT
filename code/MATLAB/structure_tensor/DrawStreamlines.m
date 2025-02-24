function DrawStreamlines(Paths,Unflipped,Flipped,SurfDist,FN,Lim)

MD = ceil(max(SurfDist(:)));
CMap = parula(101);
%VoxSize = 0.2; % 200 um
VoxSize = 0.025; % 25 um
DS = 5*VoxSize;
MaxLength = 0;

figure(FN); clf; 

subplot(1,2,1);
hold on;
for i=1:length(Unflipped)
    %fprintf('i:%d\n',i);
    LP = [flipud(Paths{Unflipped(i)}{2});Paths{Unflipped(i)}{1}];
    LPIdx = sub2ind(size(SurfDist),min(size(SurfDist,1),max(1,round(LP(:,1)))),min(size(SurfDist,2),max(1,round(LP(:,2)))),min(size(SurfDist,3),max(1,round(LP(:,3)))));
%    scatter3(LP(:,1),LP(:,2),LP(:,3),5,[0:1/(size(LP,1)-1):1]','.');
%    scatter3(LP(:,1),LP(:,2),LP(:,3),5,SurfDist(LPIdx),'.');
    if length(LPIdx)-1 > MaxLength
       MaxLength = length(LPIdx)-1 ;
    end
    XSeg = ([LP(1:end-1,1),LP(2:end,1)]-1)*VoxSize;
    YSeg = ([LP(1:end-1,2),LP(2:end,2)]-1)*VoxSize;
    ZSeg = ([LP(1:end-1,3),LP(2:end,3)]-1)*VoxSize;
    CVec = 0.5*(SurfDist(LPIdx(1:end-1))+SurfDist(LPIdx(2:end)));
    hu = plot3(XSeg',YSeg',ZSeg','-','Visible','Off');
    segColours = CMap(fix(100*(CVec./MD))+1,:);
    set(hu,{'Color'},mat2cell(segColours,ones(size(CVec)),3));
    set(hu,'Visible','on');
end
for i=1:length(Flipped)
    LP = flipud([flipud(Paths{Flipped(i)}{2});Paths{Flipped(i)}{1}]);
    LPIdx = sub2ind(size(SurfDist),min(size(SurfDist,1),max(1,round(LP(:,1)))),min(size(SurfDist,2),max(1,round(LP(:,2)))),min(size(SurfDist,3),max(1,round(LP(:,3)))));
%    scatter3(LP(:,1),LP(:,2),LP(:,3),5,[0:1/(size(LP,1)-1):1]','.');
%    scatter3(LP(:,1),LP(:,2),LP(:,3),5,SurfDist(LPIdx),'.');
    if length(LPIdx)-1 > MaxLength
        MaxLength = length(LPIdx)-1 ;
    end
    XSeg = ([LP(1:end-1,1),LP(2:end,1)]-1)*VoxSize;
    YSeg = ([LP(1:end-1,2),LP(2:end,2)]-1)*VoxSize;
    ZSeg = ([LP(1:end-1,3),LP(2:end,3)]-1)*VoxSize;
    CVec = 0.5*(SurfDist(LPIdx(1:end-1))+SurfDist(LPIdx(2:end)));
    hf = plot3(XSeg',YSeg',ZSeg','-','Visible','Off');
    segColours = CMap(fix(100*(CVec./MD))+1,:);
    set(hf,{'Color'},mat2cell(segColours,ones(size(CVec)),3));
    set(hf,'Visible','on');
end
hc=colorbar;
Delta = VoxSize*(MD/4);
set(hc,'Ticks',[0,0.25,0.5,0.75,1],'TickLabels',{'0',sprintf('%3.1f',Delta),sprintf('%3.1f',2*Delta),sprintf('%3.1f',3*Delta),sprintf('%3.1f',4*Delta)});
ylabel(hc,'Depth (mm)');
axis equal;
box on;
hold off;
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
view(3);
%view(180,0);
hold off;
if ~isempty(Lim)
  axis(Lim);
end;

% streamlines with length
subplot(1,2,2);
hold on;
for i=1:length(Unflipped)
    LP = [flipud(Paths{Unflipped(i)}{2});Paths{Unflipped(i)}{1}];
    LPIdx = sub2ind(size(SurfDist),min(size(SurfDist,1),max(1,round(LP(:,1)))),min(size(SurfDist,2),max(1,round(LP(:,2)))),min(size(SurfDist,3),max(1,round(LP(:,3)))));
%    scatter3(LP(:,1),LP(:,2),LP(:,3),5,[0:1/(size(LP,1)-1):1]','.');
%    scatter3(LP(:,1),LP(:,2),LP(:,3),5,SurfDist(LPIdx),'.');
    XSeg = ([LP(1:end-1,1),LP(2:end,1)]-1)*VoxSize;
    YSeg = ([LP(1:end-1,2),LP(2:end,2)]-1)*VoxSize;
    ZSeg = ([LP(1:end-1,3),LP(2:end,3)]-1)*VoxSize;
    CVec = [0.5:1:(length(LPIdx)-1)]';
    hu = plot3(XSeg',YSeg',ZSeg','-','Visible','Off');
    segColours = CMap(fix(100*(CVec./MaxLength))+1,:);
    set(hu,{'Color'},mat2cell(segColours,ones(size(CVec)),3));
    set(hu,'Visible','on');
end
for i=1:length(Flipped)
    LP = flipud([flipud(Paths{Flipped(i)}{2});Paths{Flipped(i)}{1}]);
    LPIdx = sub2ind(size(SurfDist),min(size(SurfDist,1),max(1,round(LP(:,1)))),min(size(SurfDist,2),max(1,round(LP(:,2)))),min(size(SurfDist,3),max(1,round(LP(:,3)))));
%    scatter3(LP(:,1),LP(:,2),LP(:,3),5,[0:1/(size(LP,1)-1):1]','.');
%    scatter3(LP(:,1),LP(:,2),LP(:,3),5,SurfDist(LPIdx),'.');
    XSeg = ([LP(1:end-1,1),LP(2:end,1)]-1)*VoxSize;
    YSeg = ([LP(1:end-1,2),LP(2:end,2)]-1)*VoxSize;
    ZSeg = ([LP(1:end-1,3),LP(2:end,3)]-1)*VoxSize;
    CVec = [0.5:1:(length(LPIdx)-1)]';
    hf = plot3(XSeg',YSeg',ZSeg','-','Visible','Off');
    segColours = CMap(fix(100*(CVec./MaxLength))+1,:);
    set(hf,{'Color'},mat2cell(segColours,ones(size(CVec)),3));
    set(hf,'Visible','on');
end
hc=colorbar;
Delta = VoxSize*(MaxLength/4);
set(hc,'Ticks',[0,0.25,0.5,0.75,1],'TickLabels',{'0',sprintf('%3.1f',Delta),sprintf('%3.1f',2*Delta),sprintf('%3.1f',3*Delta),sprintf('%3.1f',4*Delta)});
ylabel(hc,'Path Length (mm)');
axis equal;
box on;
hold off;
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
view(3);
%view(180,0);
hold off;
if ~isempty(Lim)
  axis(Lim);
end;

