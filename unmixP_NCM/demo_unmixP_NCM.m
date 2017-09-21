%% demo using sampled NCM to estimate proportion on Pavia University dataset.
% Author: Alina Zare et al. Rewriten by Sheng Zou
% Department of Electrical and Computer Engineering, University of Florida
% 09/07/2017
clear;close all;clc
load('paviaU_subimg.mat')
data = EF_reshape(paviaU_subimg);
dim = size(data,1);
M = 4; % number of endmembers
[E] = VCA(data,'Endmembers',M);
figure;plot(E);xlim([0 dim])

% cov(:,:,1)= generateSPDmatrix(dim); % randomly generate a full covariance for each endmember
% cov(:,:,2)= generateSPDmatrix(dim); % randomly generate a full covariance for each endmember
% cov(:,:,3)= generateSPDmatrix(dim); % randomly generate a full covariance for each endmember
% cov(:,:,4)= generateSPDmatrix(dim); % randomly generate a full covariance for each endmember

cov(:,:,1)= 1*eye(dim); % randomly generate a diagonal and isotropic covariance for each endmember
cov(:,:,2)= 2*eye(dim); % randomly generate a diagonal and isotropic covariance for each endmember
cov(:,:,3)= 3*eye(dim); % randomly generate a diagonal and isotropic covariance for each endmember
cov(:,:,4)= 4*eye(dim); % randomly generate a diagonal and isotropic covariance for each endmember

[Parameters] = unmixP_NCM_Parameters;
tic
[P] = unmixP_NCM(data,E,cov,Parameters);
toc

for t = 1:4
figure
imagesc(reshape(P(t,:),50,50));axis image;colorbar
end