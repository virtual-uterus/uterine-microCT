function [Wxx,Wyy,Wzz,Wxy,Wxz,Wyz] = Construct2ndDerivativeFilters(Ni,Nj,Nk,TW)

  % This function computes 2nd derivative templates of width TW (odd), and places them
  % in the lower-left-hand corner of a zero
  % array of size (Ni+4)x(Nj+4)x(Nk+4) and then reshapes that into a 1D vector.
  % In order to use the vector for an fft convolution, it is circularly 
  % shifted by -(2*(Ni+4)*(Nj+4)+2*(Ni+4)+3) and then flipped upside down
  % so that W(1) is the center of the template, W(1) is the next weight
  % moving down the original vector and W(N) is the next weight moving
  % up the original vector. Most of the W vector is zero.

  % PadWidth will be one less than half the width of the derivative template
  PadWidth = (TW-1)/2;

  % set up derivative filter - based on a generalized finite difference approximation
  [i,j,k]=ndgrid(-PadWidth:PadWidth,-PadWidth:PadWidth,-PadWidth:PadWidth);
  i = reshape(i,[prod(size(i)),1]);
  j = reshape(j,[prod(size(j)),1]);
  k = reshape(k,[prod(size(k)),1]);
  W = [i,j,k,0.5*i.^2,0.5*j.^2,0.5*k.^2,i.*j,i.*k,j.*k];
  M = pinv(W);
  M(:,((prod(size(i))-1)/2+1)) = -1*sum(M,2);

  % nxnxn filter for 2nd derivatives
  wxx = reshape(M(4,:),[TW,TW,TW]);
  wyy = reshape(M(5,:),[TW,TW,TW]);
  wzz = reshape(M(6,:),[TW,TW,TW]);
  wxy = reshape(M(7,:),[TW,TW,TW]);
  wxz = reshape(M(8,:),[TW,TW,TW]);
  wyz = reshape(M(9,:),[TW,TW,TW]);

  if TW == 3 
    wxx = padarray(wxx,[1,1,1],0,'both'); 
    wyy = padarray(wyy,[1,1,1],0,'both'); 
    wzz = padarray(wzz,[1,1,1],0,'both'); 
    wxy = padarray(wxy,[1,1,1],0,'both'); 
    wxz = padarray(wxz,[1,1,1],0,'both'); 
    wyz = padarray(wyz,[1,1,1],0,'both'); 
    TW = 5; 
    PadWidth = (TW-1)/2;
  end;

  % Circular shift offset - this is the distance to the middle 
  % of the 5x5x5 template
  CS = -(2*(Ni+2*PadWidth)*(Nj+2*PadWidth)+2*(Ni+2*PadWidth)+3);

  % effectively this is trying to build the first column of the 
  % circulant matrix (a special form of the Topelitz matrix)
  % xx derivative weights
  Wxx = zeros([Ni+2*PadWidth,Nj+2*PadWidth,Nk+2*PadWidth]);
  Wxx(1:TW,1:TW,1:TW) = wxx;
  % spread the derivative to 1D vector 
  Wxx = reshape(Wxx,[(Ni+2*PadWidth)*(Nj+2*PadWidth)*(Nk+2*PadWidth),1]);
  % shift the vector to center it at 1
  Wxx = flipud(circshift(Wxx,CS));

  % yy derivative weights
  Wyy = zeros([Ni+2*PadWidth,Nj+2*PadWidth,Nk+2*PadWidth]);
  Wyy(1:TW,1:TW,1:TW) = wyy;
  Wyy = reshape(Wyy,[(Ni+2*PadWidth)*(Nj+2*PadWidth)*(Nk+2*PadWidth),1]);
  Wyy = flipud(circshift(Wyy,CS));

  % zz derivative weights
  Wzz = zeros([Ni+2*PadWidth,Nj+2*PadWidth,Nk+2*PadWidth]);
  Wzz(1:TW,1:TW,1:TW) = wzz;
  Wzz = reshape(Wzz,[(Ni+2*PadWidth)*(Nj+2*PadWidth)*(Nk+2*PadWidth),1]);
  Wzz = flipud(circshift(Wzz,CS));

  % xy derivative weights
  Wxy = zeros([Ni+2*PadWidth,Nj+2*PadWidth,Nk+2*PadWidth]);
  Wxy(1:TW,1:TW,1:TW) = wxy;
  Wxy = reshape(Wxy,[(Ni+2*PadWidth)*(Nj+2*PadWidth)*(Nk+2*PadWidth),1]);
  Wxy = flipud(circshift(Wxy,CS));

  % xz derivative weights
  Wxz = zeros([Ni+2*PadWidth,Nj+2*PadWidth,Nk+2*PadWidth]);
  Wxz(1:TW,1:TW,1:TW) = wxz;
  Wxz = reshape(Wxz,[(Ni+2*PadWidth)*(Nj+2*PadWidth)*(Nk+2*PadWidth),1]);
  Wxz = flipud(circshift(Wxz,CS));

  % yz derivative weights
  Wyz = zeros([Ni+2*PadWidth,Nj+2*PadWidth,Nk+2*PadWidth]);
  Wyz(1:TW,1:TW,1:TW) = wyz;
  Wyz = reshape(Wyz,[(Ni+2*PadWidth)*(Nj+2*PadWidth)*(Nk+2*PadWidth),1]);
  Wyz = flipud(circshift(Wyz,CS));

return;
