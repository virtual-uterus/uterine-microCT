function [R,G,B] = GenCircSpectrum(MinVal,MaxVal)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

  RD = 0:0.166659:1;
  NIntervals = length(RD)-1;
  D = MaxVal-MinVal;
  S = MinVal;
  
  R(1,:) = [RD(1)*D+S,RD(2)*D+S,1,1];
  R(2,:) = [RD(2)*D+S,RD(3)*D+S,1,0];
  R(3,:) = [RD(3)*D+S,RD(4)*D+S,0,0];
  R(4,:) = [RD(4)*D+S,RD(5)*D+S,0,0];
  R(5,:) = [RD(5)*D+S,RD(6)*D+S,0,1];
  R(6,:) = [RD(6)*D+S,RD(7)*D+S,1,1];

  G(1,:) = [RD(1)*D+S,RD(2)*D+S,0,1];
  G(2,:) = [RD(2)*D+S,RD(3)*D+S,1,1];
  G(3,:) = [RD(3)*D+S,RD(4)*D+S,1,1];
  G(4,:) = [RD(4)*D+S,RD(5)*D+S,1,0];
  G(5,:) = [RD(5)*D+S,RD(6)*D+S,0,0];
  G(6,:) = [RD(6)*D+S,RD(7)*D+S,0,0];

  B(1,:) = [RD(1)*D+S,RD(2)*D+S,0,0];
  B(2,:) = [RD(2)*D+S,RD(3)*D+S,0,0];
  B(3,:) = [RD(3)*D+S,RD(4)*D+S,0,1];
  B(4,:) = [RD(4)*D+S,RD(5)*D+S,1,1];
  B(5,:) = [RD(5)*D+S,RD(6)*D+S,1,1];
  B(6,:) = [RD(6)*D+S,RD(7)*D+S,1,0];


end

