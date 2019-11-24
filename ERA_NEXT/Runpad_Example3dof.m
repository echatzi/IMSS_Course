%======== MATLAB Demo ====================================================
% Eleni N. Chatzi, Minas Spiridonakos, Institute of Sructural Engineering,
% ETH Zurich                                                03.06.2013
%======== Define the 3dof System (linear) ================================  
% Perform ERA on impulse (free response) data of 3dof system
clear all;close all;clc;

% Mass Matrix
M=.001*[1 0 0;0 1 0;0 0 1];

% Stiffness Matrix
K=[4 -2 0;-2 4 -2;0 -2 2];

% Damping ratio  (\zeta_i)
xi=0.002;

%======== Natural Frequencies, Modes ====================== 
% Solve the eigenvalue problem
[V,D]=eig(K,M); %V:eigenvectors, D:eigenvalues(=w^2)
% Note: The Eigenvectors V are already Mass Normalized in MATLAB

% Natural Frequencies
w=[sqrt(D(1,1)) sqrt(D(2,2)) sqrt(D(3,3))];

%======== Define Raleigh damping ===========================
% alpha and beta
beta=2*xi/(w(1)+w(2));
alpha=2*xi*w(1)-beta*w(1)^2;

% Damping Matrix
C=alpha*M+beta*K;

% Discrete Time Domain Info
fs=100;                  %Sampling Frequency (has to be above 2*max freq expected in response signal)
dt=1/fs;                %Sampling interval
Ttot=20;                %Total analysis time in seconds
time=[0:dt:Ttot];       %Time vector for this example
N=length(time);         %Number of points
 
%======== Define the type of Excitation ===========================
inptype = 'known';      %Options 'imp'/'WN'/'known'

switch inptype
%%%%%%%%%%%%%%% A. UNIT IMPULSE %%%%%%%%%%%%%%%
case 'imp'
f=zeros(N,1);
f(1)=1;
ref = [];
%%%%%%%%%%%%%%% B. White noise excitation (ambient-unmeasured) %%%%%%%%%%%%%%%
case 'WN'      %%Choose this option if input is white noise and you desire otput only id, 
% exploiting the fact that the cross correlation of response signals to any reference one satisfy the free response equation
f=1*randn(N,1); 
ref=2;      %Determine the reference channel, i.e. the channel with respect to which you will obtain the cross correlation function%%%%%%%%%%%%%%% C. Measured random excitation (could also be earthquake) %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% C. Aribtrary but known excitation %%%%%%%%%%%%%%%
case 'known'
load ElCentro.mat
f=1000*ElCentro(1:N,2);
ref = [];
end

%% ======== Continuous State Space Form ===========================
Ac = [zeros(3) eye(3);-inv(M)*K -inv(M)*C];
% print modal quantities for checking
[Wn,zeta] = damp(Ac);

display('Frequencies of the real system (Hz)');
Wn([6 4 2])/2/pi

Bc = [zeros(3,1);diag(inv(M))];
Cc = [-inv(M)*K -inv(M)*C];            %assuming we measure accelerations
Dc = zeros(3,1);

sys0=ss(Ac,Bc,Cc,Dc);
[Y,T,X]=lsim(sys0,f,time,zeros(6,1));
output=Y;                     %Here you choose which measurement channels you would like to use

% Plot the spectrum of the 3rd DOF response 
figure
pwelch(Y(:,3),[],[],[],fs);
pause
%% ======== call the ERA ===========================
nch=size(Y,2);
Nfft = 2^(nextpow2(N)-1);
 
%% Call The ERA
% Specify the number of modes you are interested in identifying
% Run for increasing # ndof 
[Pxx,fxx] = pwelch(Y(:,3),[],[],[],fs);
plot(fxx,log(Pxx)/4);hold on

range = [2:6];
for ndof = range
 %Define the size of the Hankel matrix - crucial for convergence, depends on data
 order=4*ndof;      %Recommended value = 4*number of modes

[freq,err] = ERA_NEXT_3dof(output,f,inptype,nch,ref,ndof,fs,Nfft,order);
if err~=0
    break
end
scatter(freq,ndof*ones(ndof,1),'filled');hold on;grid on
end
xlabel('Frequency (Hz)','fontweight','bold','fontsize',14)
ylabel('Number of Modes','fontweight','bold','fontsize',14)
title ('Stabilization Chart','fontweight','bold','fontsize',14)
pause

%% Select the Optimal ndof and run the ERA for this configuration
ndof=input('Select Optimum Order: ndof= ');
order=4*ndof;  %20 for imp/'known'
[freq_rel2,err,phi,ao,bo,co] = ERA_NEXT_3dof(output,f,inptype,nch,ref,ndof,fs,Nfft,order);

%% Plot identified versus true (reference) modeshapes
figure
for i=1:ndof
subplot(1,ndof,i)
[mval, ind]=max(abs(phi(:,i)));
p1=plot([0:ndof],[0 phi(:,i)'*sign(phi(ind,i))/mval]);hold on;grid on
[mval, ind]=max(abs(V(:,i)));
p2=plot([0:3],[0 V(:,i)'*sign(V(ind,i))/mval],'r');
legend([p1 p2],'calculated','analytical')
title (['Mode ',num2str(i)]);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reproduce the system response for a random input, using the ERA estimated matrices
% This is only possible for identification under availability of impulse response, otherwise the bo matrix is not appropriately identified.
switch inptype
case 'imp'
fnew=randn(N,1);

[Y1,T1,X1]=lsim(sys0,fnew,time,zeros(2*ndof,1));
sys=ss(ao,bo,co,zeros(ndof,1),dt);
[Y2,T2,X2]=lsim(sys,fnew,time,zeros(2*ndof,1));

figure
for i=1:ndof
subplot(ndof,1,i)
plot(Y2(:,i));hold on;
plot(Y1(:,i),'r--');grid on
legend('ERA estimate', 'true response')
title(['Response at Measured dof ',num2str(i)])
end
end



