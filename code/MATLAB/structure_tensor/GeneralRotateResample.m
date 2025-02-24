function [S,Xs,Ys,Zs] = GeneralRotateResample(I,location,DMin,DMax,f,s,n)

%function [S,Xs,Ys,Zs] = GeneralRotateResample(I,location,D,f,s,n)
  % This function resamples mask image I at location
  % to a block D1xD2xD3 voxels oriented so that f, s and n
  % are orthogonal to the block.

  % Array of offsets (force all Ds to be odd)
%   if ~mod(D(1),2) D(1) = D(1)-1; end
%   if ~mod(D(2),2) D(2) = D(2)-1; end
%   if ~mod(D(3),2) D(3) = D(3)-1; end
%   d = floor((D-1)/2); % D is odd
%   [F,S,N]=ndgrid(-d(1):d(1),-d(2):d(2),-d(3):d(3)); 
%  [F,S,N]=meshgrid(-d(1):d(1),-d(2):d(2),-d(3):d(3)); 
 % [F,N,S]=ndgrid(-d(1):d(1),d(2):-1:-d(2),-d(3):d(3)); % Locate interpolated image with origin in upper left

  [F,S,N]=ndgrid(DMin(1):DMax(1),DMin(2):DMax(2),DMin(3):DMax(3)); 

  % Rotated image sample points
  Xs = location(1)+F*f(1)+S*s(1)+N*n(1);
  Ys = location(2)+F*f(2)+S*s(2)+N*n(2);
  Zs = location(3)+F*f(3)+S*s(3)+N*n(3);

  % Set ranges on I
  Rmin = max(floor([min(min(min(Xs))),min(min(min(Ys))),min(min(min(Zs)))])-[1,1,1],[1,1,1]); 
  Rmax = min(floor([max(max(max(Xs))),max(max(max(Ys))),max(max(max(Zs)))])+[1,1,1],size(I)); 

  % Resample image  
%  S = iminterpn(single(I(Rmin(1):Rmax(1),Rmin(2):Rmax(2),Rmin(3):Rmax(3))),Xs-Rmin(1)+1,Ys-Rmin(2)+1,Zs-Rmin(3)+1,'linear',0);
%  S = iminterpn(I(Rmin(1):Rmax(1),Rmin(2):Rmax(2),Rmin(3):Rmax(3)),Xs-Rmin(1)+1,Ys-Rmin(2)+1,Zs-Rmin(3)+1,'nearest',0);
%  S = iminterpn(I,Xs-Rmin(1)+1,Ys-Rmin(2)+1,Zs-Rmin(3)+1,'linear',0);
  
  if strcmp(class(I),'logical')
      S = iminterpn(I,Xs,Ys,Zs,'nearest',0);
  else
      S = iminterpn(I,Xs,Ys,Zs,'linear',0);
  end

return;
