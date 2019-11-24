function [freq_rel2,err,phi,ao,bo,co] = ERA_NEXT_3dof(output,f,inptype,nch,ref,ndof,fs,Nfft,order)

%=============================================================================================
% This is the main file for calling the ERA function (ERA.m)

switch inptype
    case 'imp'
     dt=1/fs;   
     YY=output;      %This is already impulse response
    
    case 'WN' 
    maxlag=length(f);
    for ii=1:nch
        Rxy(:,ii) = xcorr(output(:,ii),output(:,ref),maxlag,'coeff');
        [csdxy(:,ii), Fxy]= cpsd(output(:,ref),output(:,ii),[],[],Nfft,fs); %CSD of the output and input1.
  
    % Mirroring with complex conjugate
    df = Fxy(2)-Fxy(1);
    FF = 0:df:Nfft*df;
    CSD(:,ii)=[csdxy(:,ii); conj(flipud(csdxy(1:end-1,ii)))]; % Be sure not to overlap the data.
    dt = 1/df/Nfft;
    TT = 0: dt : (Nfft-1)*dt;
    YY(:,ii) = 1/dt * real(ifft(CSD(:,ii), Nfft));
    
        % Truncating to positive lags
        RY(:,ii) = Rxy(maxlag+2:end,ii);
 
    end
   
    otherwise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Impulse Response Function for ERA (NExT method) in case of specified
%(known) input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalize the input force and the signals 

input=f;
norm=1;
input_norm=input/norm;
output_norm=output/norm;

    [csdxx,Fxx] = cpsd(input_norm,input_norm,[],[],Nfft,fs); %PSD of the input1 signal.
    for ii = 1:nch
        [csdxy(:,ii), Fxy]= cpsd(input_norm,output_norm(:,ii),[],[],Nfft,fs); %CSD of the output and input1.
    end

    % Obtain the transfer function from H = Sxy/Sxx.
    for ii = 1:nch
        Hxy(:,ii) = csdxy(:,ii)./csdxx;
    end

%     Mirroring with complex conjugate
    for ii = 1:nch
        Htemp(:,ii)=[Hxy(:,ii); conj(flipud(Hxy(1:end-1,ii)))]; % Be sure not to overlap the data.
    end
    df = Fxy(2)-Fxy(1);
    FF = 0:df:Nfft*df;

%     Impulse response function
    dt = 1/df/Nfft;
    TT = 0: dt : (Nfft-1)*dt;
    for ii = 1:nch
        YY(:,ii) = 1/dt * real(ifft(Htemp(:,ii), Nfft));
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
[y_ident_full,yx,freq_rel,damp_rel,modal_shapes,freq_rel2,ao,bo,co] = ...
    ERA([YY],dt,1,1,size(YY,1),2*ndof,order);

err=0;
catch ME
    if ME~=0
    err=1;
    freq_rel2=0;
    ao=0;bo=0;co=0;phi=0;
    return
    end
end

% Print/Store Outputs
[freq_rel2 ind]=sort(freq_rel2);           %sort frequencies starting from lowest one

display('Identified Frequencies (Hz)');
freq_rel2 = freq_rel2(1:2:end)
phi=real(modal_shapes(:,ind(1:2:end)));  %sort the eigenvectors correspondingly, and only keep them once
