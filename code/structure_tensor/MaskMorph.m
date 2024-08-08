 function Mask = MaskMorph(Mask1,Mask2,L1,L2,L)

  % This function morphs between two logical masks.
  % Input:
  %       Mask1 - first mask, located at L1
  %       Mask2 - second mark, location at L2
  %       L - evaluation location between L1 and L2
  % Output:
  %       Mask - a morphed mask
  %
  % Example:
  %   Mask = MaskMorph(Mask1,Mask2,0,2,1);
  %tstart = tic;

  alpha = (L-L1)/(L2-L1);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % Areas in Mask2 and are not in Mask1
  M2M1 = false(size(Mask1));
  M2M1(setdiff(find(Mask2),find(Mask1))) = 1;
  
  % Initial distance measure based on Mask1 
  Dist1 = single(M2M1).*bwdist(Mask1); 

  % Areas that are not adjacent to Mask1 will have a minimum distance
  % greater than 1
  stats = regionprops(M2M1,'PixelIdxList','Centroid');
  MinDist = arrayfun(@(x)min(Dist1(x.PixelIdxList)),stats);
  
  % When applying to full images a lot of centroids are activated for some
  % reason. Just ignore the non-adjacent feature - not important and will
  % save operations.
  % Areas that are not adjacent should have distances growing from the
  % centroid.
%  IdxNonAdjacent = find(MinDist > 1);
%  if ~isempty(IdxNonAdjacent),
%    Centroid = round(cat(1,stats(IdxNonAdjacent).Centroid));
%    Mask1(Centroid(:,2),Centroid(:,1)) = 1;
%    Dist1 = single(M2M1).*bwdist(Mask1);
%  end;
%  if ~isempty(IdxNonAdjacent),
%    MaskWork = Mask1;
%    for i=1:length(IdxNonAdjacent),
%      Centroid = round(stats(IdxNonAdjacent(i)).Centroid);
%      MaskWork(Centroid(2),Centroid(1)) = 1;
%    end;
%    Dist1 = bwdist(MaskWork);
%  end;
%  Dist1 = single(M2M1).*Dist1;
  
  % Normalize the distance of each area
  for i=1:length(stats),
      Dist1(stats(i).PixelIdxList) = Dist1(stats(i).PixelIdxList)/...
            max(Dist1(stats(i).PixelIdxList));
  end;

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Repeat same process with Dist2 etc - that will be used to grow Mask2 to
  % Mask1
  
  % Areas in Mask1 and are not in Mask2
  M1M2 = false(size(Mask2));
  M1M2(setdiff(find(Mask1),find(Mask2))) = 1;

  % Initial distance measure based on Mask2 
  Dist2 = single(M1M2).*bwdist(Mask2); 

  % Areas that are not adjacent to Mask2 will have a minimum distance
  % greater than 1
  stats = regionprops(M1M2,'PixelIdxList','Centroid');
  MinDist = arrayfun(@(x)min(Dist2(x.PixelIdxList)),stats);
  
  % Areas that are not adjacent should have distances growing from the
  % centroid.
%  IdxNonAdjacent = find(MinDist > 1);
%  if ~isempty(IdxNonAdjacent),
%    Centroid = round(cat(1,stats(IdxNonAdjacent).Centroid));
%    Mask2(Centroid(:,2),Centroid(:,1)) = 1;
%    Dist2 = single(M1M2).*bwdist(Mask2);
%  end;
%  if ~isempty(IdxNonAdjacent),
%    MaskWork = Mask2;
%    for i=1:length(IdxNonAdjacent),
%      Centroid = round(stats(IdxNonAdjacent(i)).Centroid);
%      MaskWork(Centroid(2),Centroid(1)) = 1;
%    end;
%    Dist2 = bwdist(MaskWork);
%  end;
%  Dist2 = single(M1M2).*Dist2;
  
  % Normalize the distance of each area
  for i=1:length(stats),
    Dist2(stats(i).PixelIdxList) = Dist2(stats(i).PixelIdxList)/...
            max(Dist2(stats(i).PixelIdxList));
  end;

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Find intermediate mask by advancing Mask1 by a proportion (alpha) of 
  % Dist1 (OR) and removing a proportion (1-alpha) of Dist2 (XOR)
  Mask = xor((Mask1 | (Dist1 > 0 & Dist1 <= alpha)),(Dist2 > (1-alpha)));
 
  %telapsed = toc(tstart);

  %fprintf('Time: %f \n',telapsed);
return;
