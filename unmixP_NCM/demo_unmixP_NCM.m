%% demo using sampled NCM to estimate proportion on Pavia University dataset.
% Author: Alina Zare et al. Rewriten by Sheng Zou
% Department of Electrical and Computer Engineering, University of Florida
% 09/07/2017
clear;close all;clc
load('paviaU_subimg.mat') % load a subset of Pavia University hyperspectral image, the size is 50-by-50-by-103
data = EF_reshape(paviaU_subimg);% reshape the 3D HSI data to 2D (d-by-(M*N))
dim = size(data,1); % get the dimensionality d of the HSI
M = 4; % preset number of endmembers
[E] = VCA(data,'Endmembers',M);% Run VCA (Vertex component analysis) to estimate the endmember signatures for the HSI data
figure;plot(E);xlim([0 dim]) % show the estimated endmember signatures

% cov(:,:,1)= generateSPDmatrix(dim); % randomly generate a full covariance for each endmember
% cov(:,:,2)= generateSPDmatrix(dim); % randomly generate a full covariance for each endmember
% cov(:,:,3)= generateSPDmatrix(dim); % randomly generate a full covariance for each endmember
% cov(:,:,4)= generateSPDmatrix(dim); % randomly generate a full covariance for each endmember

cov(:,:,1)= 1*eye(dim); % hand pick a diagonal and isotropic covariance for the endmember
cov(:,:,2)= 2*eye(dim); % hand pick a diagonal and isotropic covariance for the endmember
cov(:,:,3)= 3*eye(dim); % hand pick a diagonal and isotropic covariance for the endmember
cov(:,:,4)= 4*eye(dim); % hand pick a diagonal and isotropic covariance for the endmember

[Parameters] = unmixP_NCM_Parameters; % load the parameters associated with sampled NCM algorithm
tic
[P] = unmixP_NCM(data,E,cov,Parameters); % run the main function of sampled NCM to estimate the proportions 
toc

for t = 1:4
figure
imagesc(reshape(P(t,:),50,50));axis image;colorbar % show each proportion map
end