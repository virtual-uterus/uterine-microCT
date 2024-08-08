function [ST] = ComputeStructureTensor(I,DTW)
    % ComputeStructureTensor
    % This function finds the derivatives of 3D data in I and determines the
    % structure tensor components
    %
    % Inputs: I - data; real (single or double), normalised 0-1
    %         DTW - derivative template width
    % Outputs: ST - 6 entry cell array with arrays of ST values

    % Sizes and steps
    [Ni,Nj,Nk]=size(I);
    PD = 2; % the templates are always 5x5x5
    %PD = floor(DTW/2);

    % Pad the image around each edge - reflect image in padding
    %fprintf('... padding array ...\n');
    IPad = padarray(I,[PD,PD,PD],'symmetric'); clear I;
    IPad = reshape(IPad,[(Ni+2*PD)*(Nj+2*PD)*(Nk+2*PD),1]);

    % Reshape padded permuted equalized image to a 1D array
    %fprintf('... computing FFT of image ...\n');
    IPadf = fft(IPad); clear IPad;

    % Farid and Simoncelli derivative filter
    if DTW == 5
        % spread weights
        p = [0.037659  0.249153  0.426375  0.249153  0.037659];
        % derivative weights
        d1 =[0.109604  0.276691  0.000000 -0.276691 -0.109604];
    elseif DTW == 3
        % spread weights
        p = [0.0  0.229879  0.540242  0.229879  0.0];
        % derivative weights
        d1 =[0.0  0.425287  0.000000 -0.425287 -0.0];
    else
        fprintf(' %%% This derivative template width is not implemented! %%%\n');
        return;
    end

    % complete nxnxn filter for derivatives in z
    wz = repmat(p'*p,[1,1,5]);
    wz(:,:,1) = wz(:,:,1)*d1(5);  
    wz(:,:,2) = wz(:,:,2)*d1(4);  
    wz(:,:,3) = wz(:,:,3)*d1(3);  
    wz(:,:,4) = wz(:,:,4)*d1(2);  
    wz(:,:,5) = wz(:,:,5)*d1(1);
    wz = cast(wz,class(IPadf));

    % set derivative templates on other directions by permuting the dz
    % template
    wx = permute(wz,[3,1,2]);
    wy = permute(wz,[1,3,2]);

    % Circular shift offset - this is the distance to the middle 
    % of the 5x5x5 template
    CS = -(2*(Ni+2*PD)*(Nj+2*PD)+2*(Ni+2*PD)+3);
    
    % Compute fft of 1st derivative weights
    %fprintf('... computing FFTs of 1st derivative weights ...\n');
    %fprintf('... performing FFT convolution for 1st derivatives ...\n');

    % effectively this is trying to build the first column of the 
    % circulant matrix (a special form of the Topelitz matrix)
    % x derivative weights
    Wx = cast(false([Ni+2*PD,Nj+2*PD,Nk+2*PD]),class(wx));
    Wx(1:5,1:5,1:5) = wx;
    % spread the derivative to 1D vector 
    Wx = reshape(Wx,[(Ni+2*PD)*(Nj+2*PD)*(Nk+2*PD),1]);
    % shift the vector to center it at 1
    Wx = flipud(circshift(Wx,CS));
    % FFT of derivative weight
    Wxf = fft(Wx); clear Wx;
    % Derivative with inverse FFT
    Di = ifft(IPadf .* Wxf); clear Wxf;% FFT convolution
    Di = reshape(Di,[Ni+2*PD,Nj+2*PD,Nk+2*PD]);
    Di = Di(1+PD:Ni+PD,1+PD:Nj+PD,1+PD:Nk+PD);

    % y derivative weights
    Wy = cast(zeros([Ni+2*PD,Nj+2*PD,Nk+2*PD]),class(wy));
    Wy(1:5,1:5,1:5) = wy;
    Wy = reshape(Wy,[(Ni+2*PD)*(Nj+2*PD)*(Nk+2*PD),1]);
    Wy = flipud(circshift(Wy,CS));
    % FFT of derivative weight
    Wyf = single(fft(Wy)); clear Wy;
    % Derivative with inverse FFT
    Dj = ifft(IPadf .* Wyf); clear Wyf;% FFT convolution
    Dj = reshape(Dj,[Ni+2*PD,Nj+2*PD,Nk+2*PD]);
    Dj = Dj(1+PD:Ni+PD,1+PD:Nj+PD,1+PD:Nk+PD);

    % z derivative weights
    Wz = cast(zeros([Ni+2*PD,Nj+2*PD,Nk+2*PD]),class(wz));
    Wz(1:5,1:5,1:5) = wz;
    Wz = reshape(Wz,[(Ni+2*PD)*(Nj+2*PD)*(Nk+2*PD),1]);
    Wz = flipud(circshift(Wz,CS));
    % FFT of derivative weight
    Wzf = single(fft(Wz)); clear Wz;
    % Derivative with inverse FFT
    Dk = ifft(IPadf .* Wzf); clear Wzf;% FFT convolution
    Dk = reshape(Dk,[Ni+2*PD,Nj+2*PD,Nk+2*PD]);
    Dk = Dk(1+PD:Ni+PD,1+PD:Nj+PD,1+PD:Nk+PD);
    
    % Structure tensor cell array
    ST = cell(6,1);
    ST{1} = Di.*Di; %ST{1} = reshape(ST{1},[Ni,Nj,Nk]); 
    ST{2} = Di.*Dj; %ST{2} = reshape(ST{2},[Ni,Nj,Nk]);
    ST{3} = Di.*Dk; %ST{3} = reshape(ST{3},[Ni,Nj,Nk]); clear Di; 
    ST{4} = Dj.*Dj; %ST{4} = reshape(ST{4},[Ni,Nj,Nk]);
    ST{5} = Dj.*Dk; %ST{5} = reshape(ST{5},[Ni,Nj,Nk]); clear Dj; 
    ST{6} = Dk.*Dk; %ST{6} = reshape(ST{6},[Ni,Nj,Nk]); clear Dk; 
end

