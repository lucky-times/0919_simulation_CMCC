% clear all;close all;
% load c:\temp\MIMO.mat

nTx=NumberOfAntennas_NodeB; % Number of Tx Antennas, M
nRx=NumberOfAntennas_UE; % Number of Rx  Antennas, N
lH=length(H); % Length of Channel
type=0; % 0=complex
% R = Correlation Matrix
% H = Channel 
% H = generate_H(nTx, nRx, R, lH, type, PDP_linear); % Eikös meillä ole jo H!!!???
% Hu(:,:,:)=H(:,:,1,:); % shape H
% Pw=water_fill(Plinear,Le,ee)
% produces matrices of eigenvalues (D) and eigenvectors (yyy) of matrix A
[yyy D]=eig(R); % 
clear yyy;
[yyy Pw]=eig(R); %
Pw2=Pw(1:nTx,1:nRx);
Pw=Pw2;
clear yyy;
% returns a vector of the eigenvalues of matrix R

V=zeros(nTx,nTx,lH);
for i=1:lH,
        [V(:,:,i) xxx]=eig(H(:,:,i)); 
        clear xxx;
end

VH=zeros(size(V));
for j=1:lH,
    VH(:,:,j)=conj(V(:,:,j));
end
 
% size(V(:,:,1))
% size(Pw)
% size(VH(:,:,1))

for k=1:lH,
    Sw(:,:,k)=V(:,:,k)*Pw*VH(:,:,k);
end


Plinear=2; % Total available power for transmission (W)
sigma_N_2=1.05; % Noise variance
e=ones(nTx,nRx);

sum_inst_cap = capacityMM(nTx, nRx, H, e, V, Plinear, sigma_N_2)