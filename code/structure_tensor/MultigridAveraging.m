function [S,SI,SJ,SK]=MultigridAveraging(D,DI,DJ,DK,TW)

  % This function takes arrays of valid data D and uses a level 4 binomial 
  % filter to smooth to the next level of coarseness.

  % Level 4 binomial filter
  if TW == 5,
    Wm = [1,4,6,4,1]/16; 
  elseif TW == 3,
    Wm = [0,1,2,1,0]/4;
  else
    fprintf(' %%% This smoothing template width is not implemented! %%%\n');
  end;
  PL = 2; % padding level for filter

  % Double step ranges and offset for smoothing
  [Ni,Nj,Nk]=size(D);
  I = [1+1*PL:2:Ni-1*PL]; Nmi = length(I); MaxI = max(I)+2;
  J = [1+1*PL:2:Nj-1*PL]; Nmj = length(J); MaxJ = max(J)+2;
  K = [1+1*PL:2:Nk-1*PL]; Nmk = length(K); MaxK = max(K)+2;
  
  % Smooth in i
  fprintf(' ... Smoothing in i.');
  t0 = clock;
  I1 = [3,I-2,MaxI-2];
  I2 = [2,I-1,MaxI-1];
  I3 = [1,I,MaxI];
  I4 = [2,I+1,MaxI-1];
  I5 = [3,I+2,MaxI-2];
  S = Wm(1)*D(I1,:,:)+Wm(2)*D(I2,:,:)+Wm(3)*D(I3,:,:)+Wm(4)*D(I4,:,:)+Wm(5)*D(I5,:,:);
  SI = DI(I3,:,:); SJ = DJ(I3,:,:); SK = DK(I3,:,:);
  t1 = etime(clock,t0); fprintf(' Time: %0.2f sec\n',t1);

  % Smooth in j
  fprintf(' ... Smoothing in j.');
  t0 = clock;
  J1 = [3,J-2,MaxJ-2];
  J2 = [2,J-1,MaxJ-1];
  J3 = [1,J,MaxJ];
  J4 = [2,J+1,MaxJ-1];
  J5 = [3,J+2,MaxJ-2];
  S = Wm(1)*S(:,J1,:)+Wm(2)*S(:,J2,:)+Wm(3)*S(:,J3,:)+Wm(4)*S(:,J4,:)+Wm(5)*S(:,J5,:);
  SI = SI(:,J3,:); SJ = SJ(:,J3,:); SK = SK(:,J3,:);
  t1 = etime(clock,t0); fprintf(' Time: %0.2f sec\n',t1);

  % Smooth in k
  fprintf(' ... Smoothing in k.');
  t0 = clock;
  K1 = [3,K-2,MaxK-2];
  K2 = [2,K-1,MaxK-1];
  K3 = [1,K,MaxK];
  K4 = [2,K+1,MaxK-1];
  K5 = [3,K+2,MaxK-2];
  S = Wm(1)*S(:,:,K1)+Wm(2)*S(:,:,K2)+Wm(3)*S(:,:,K3)+Wm(4)*S(:,:,K4)+Wm(5)*S(:,:,K5);
  SI = SI(:,:,K3); SJ = SJ(:,:,K3); SK = SK(:,:,K3);
  t1 = etime(clock,t0); fprintf(' Time: %0.2f sec\n',t1);
  
return;
