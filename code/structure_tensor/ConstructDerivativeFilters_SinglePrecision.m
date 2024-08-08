function [Wx,Wy,Wz] = ConstructDerivativeFilters_SinglePrecision(Ni,Nj,Nk,TW)

  % This function takes the "optimal" nxnxn (n=3 or 5) derivative templates as 
  % defined by Farid & Simoncelli, places it in the lower-left-hand corner of a zero
  % array of size (Ni+4)x(Nj+4)x(Nk+4) and then reshapes that into a 1D vector.
  % In order to use the vector for an fft convolution, it is circularly 
  % shifted by -(2*(Ni+4)*(Nj+4)+2*(Ni+4)+3) and then flipped upside down
  % so that W(1) is the center of the template, W(1) is the next weight
  % moving down the original vector and W(N) is the next weight moving
  % up the original vector. Most of the W vector is zero.

  % PadWidth will be one less than half the width of the derivative template
  PadWidth = 2;

  % Farid and Simoncelli derivative filter
  if TW == 5,
    % spread weights
    p = single([0.037659  0.249153  0.426375  0.249153  0.037659]);
    % derivative weights
    d1 =single([0.109604  0.276691  0.000000 -0.276691 -0.109604]);
  elseif TW == 3,
    % spread weights
    p = single([0.0  0.229879  0.540242  0.229879  0.0]);
    % derivative weights
    d1 = single([0.0  0.425287  0.000000 -0.425287 -0.0]);
  else
    fprintf(' %%% This derivative template width is not implemented! %%%\n');
  end;

  % complete nxnxn filter for derivatives in z
  wz = repmat(p'*p,[1,1,5]);
  wz(:,:,1) = wz(:,:,1)*d1(5);  
  wz(:,:,2) = wz(:,:,2)*d1(4);  
  wz(:,:,3) = wz(:,:,3)*d1(3);  
  wz(:,:,4) = wz(:,:,4)*d1(2);  
  wz(:,:,5) = wz(:,:,5)*d1(1);

  % set derivative templates on other directions by permuting the dz
  % template
  wx = permute(wz,[3,1,2]);
  wy = permute(wz,[1,3,2]);

  % Circular shift offset - this is the distance to the middle 
  % of the 5x5x5 template
  CS = -(2*(Ni+2*PadWidth)*(Nj+2*PadWidth)+2*(Ni+2*PadWidth)+3);

  % effectively this is trying to build the first column of the 
  % circulant matrix (a special form of the Topelitz matrix)
  % x derivative weights
  Wx = single(zeros([Ni+2*PadWidth,Nj+2*PadWidth,Nk+2*PadWidth]));
  Wx(1:5,1:5,1:5) = wx;
  % spread the derivative to 1D vector 
  Wx = reshape(Wx,[(Ni+2*PadWidth)*(Nj+2*PadWidth)*(Nk+2*PadWidth),1]);
  % shift the vector to center it at 1
  Wx = flipud(circshift(Wx,CS));

  % y derivative weights
  Wy = single(zeros([Ni+2*PadWidth,Nj+2*PadWidth,Nk+2*PadWidth]));
  Wy(1:5,1:5,1:5) = wy;
  Wy = reshape(Wy,[(Ni+2*PadWidth)*(Nj+2*PadWidth)*(Nk+2*PadWidth),1]);
  Wy = flipud(circshift(Wy,CS));

  % z derivative weights
  Wz = single(zeros([Ni+2*PadWidth,Nj+2*PadWidth,Nk+2*PadWidth]));
  Wz(1:5,1:5,1:5) = wz;
  Wz = reshape(Wz,[(Ni+2*PadWidth)*(Nj+2*PadWidth)*(Nk+2*PadWidth),1]);
  Wz = flipud(circshift(Wz,CS));

return;
