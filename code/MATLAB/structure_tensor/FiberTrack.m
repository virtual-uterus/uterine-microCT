function [Paths] = FiberTrack(StartLoc,DS,I,J,K,Fd2Xs,FdXYs,FdXZs,Fd2Ys, ...
    FdYZs,Fd2Zs,Mask,NM,FiberIndex,TrackLength, centre_points, nb_used_slices)

% This function interpolates the structure tensor field at a sequence of
% points along a fiber track

% FiberIndex = 1 if smallest eigenvalue direction is used otherwise it can
% be 2.

% Track 
Paths = cell(4,1);
Length = zeros(2,1);
W = [1,-1];
for i=1:2
  Length(i) = 0;
  CurrentPoint = StartLoc;
  OldFiber = [];
  Paths{i} = [CurrentPoint];
  %while Length < MaxLength

  TP = [min(NM(1),max(1,round(CurrentPoint(1)))),min(NM(2),max(1,round(CurrentPoint(2)))),min(NM(3),max(1,round(CurrentPoint(3))))];
  Test = Mask(sub2ind(NM,TP(1),TP(2),TP(3)));

  NSteps = 0;
  
  while Test & NSteps < TrackLength
      S = [Fd2Xs(CurrentPoint(1),CurrentPoint(2),CurrentPoint(3)),...
           FdXYs(CurrentPoint(1),CurrentPoint(2),CurrentPoint(3)),...
           FdXZs(CurrentPoint(1),CurrentPoint(2),CurrentPoint(3));...
           FdXYs(CurrentPoint(1),CurrentPoint(2),CurrentPoint(3)),...
           Fd2Ys(CurrentPoint(1),CurrentPoint(2),CurrentPoint(3)),...
           FdYZs(CurrentPoint(1),CurrentPoint(2),CurrentPoint(3));...
           FdXZs(CurrentPoint(1),CurrentPoint(2),CurrentPoint(3)),...
           FdYZs(CurrentPoint(1),CurrentPoint(2),CurrentPoint(3)),...
           Fd2Zs(CurrentPoint(1),CurrentPoint(2),CurrentPoint(3))];
      [V,D] = eig(S); % evect/eval in largest to smallest
      [y,idx]=sort(diag(D));
      L1 = D(idx(3),idx(3));
      L2 = D(idx(2),idx(2));
      L3 = D(idx(1),idx(1));
      NewNormal = V(:,idx(3))';
      NewSheet = V(:,idx(2))';
      NewFiber = V(:,idx(FiberIndex))';
    
      angle = ComputeFibreAngle(NewFiber, centre_points, ...
          CurrentPoint(1), CurrentPoint(3), nb_used_slices);

      if angle > 90
          angle = 180 - angle;
      end
      
      NewFiber = W(i)*NewFiber;
      NewSheet = W(i)*NewSheet;

      % Ensure that path continues moving forward in continuous direction
      if ~isempty(OldFiber)
          if dot(NewFiber,OldFiber)< 0
              NewFiber = -1*NewFiber;
              NewSheet = -1*NewSheet;
          end
      end

%       Paths{i+2} = [Paths{i+2};[L3,L2,L1]];
      
      % Find predicted/test point
      TestPoint = CurrentPoint + DS*NewFiber;
      
      % Use test point if within image range and within mask
      if norm(min(max(TestPoint,[1,1,1]),NM)-TestPoint) < 1e-6
         
          RTP = round(TestPoint);
          if Mask(sub2ind(NM,RTP(1),RTP(2),RTP(3)))
              CurrentPoint = TestPoint;
              OldFiber = NewFiber;

              % Calculate fibre angle
              angle = ComputeFibreAngle(NewFiber, centre_points, ...
                  CurrentPoint(1), CurrentPoint(3), nb_used_slices);

              if angle > 90
                  angle = 180 - angle;
              end
              Paths{i} = [Paths{i};CurrentPoint];
              Length(i) = Length(i)+DS;
              Paths{i+2} = [Paths{i+2};[L3,L2,L1, angle]];
              
              Test = 1;
          else
              Test = 0;
          end
      else
          Test = 0;
      end
      
%       DeltaS = min([(min(max(TestPoint,[1,1,1]),NM)-CurrentPoint)./NewFiber,DS]);
%       
%       CurrentPoint = CurrentPoint+DeltaS*NewFiber;
%       OldFiber = NewFiber;
%     
%       Paths{i} = [Paths{i};CurrentPoint];
%       Length(i) = Length(i)+DS;
% 
%       TP = [min(NM(1),max(1,round(CurrentPoint(1)))),min(NM(2),max(1,round(CurrentPoint(2)))),min(NM(3),max(1,round(CurrentPoint(3))))];
%       Test = Mask(sub2ind(NM,TP(1),TP(2),TP(3)));
%       
%       if DeltaS < DS Test=0; end
      
      NSteps = NSteps + 1;

  end
  Paths{i+2} = [Paths{i+2};[L3,L2,L1, angle]]; % repeat final eigenvalue content

end
