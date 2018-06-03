function knn_out = knn_classifier(hsi_img,train_data,mask,k)
%function knn_out = knn_classifier(hsi_img,train_data,mask,k)
%
%  a simple K nearest nieghbors classifier 
%
% Inputs:
%   hsi_img: hyperspectral data cube (n_rows x n_cols x n_bands)
%   train_data: structure containing training data 
%       train_data(i).Spectra: matrix containing training data from class i
%   mask: binary image indicating where to apply classifier
%   k:  number of neighbors to use during classification
%
% 10/31/2012 - Taylor C. Glenn
% 05/12/2018 - Edited by Alina Zare


if ~exist('mask','var'); mask = []; end
knn_out = img_det(@knn_cfr,hsi_img,train_data,mask,k);

end

function knn_out = knn_cfr(hsi_data,train_data,K)

% concatenate the training data
train = [train_data.Spectra];
n_train = size(train,2);

labels = zeros(n_train,1);
n_class = numel(train_data);
last = 0;
for i=1:n_class
    nt = size(train_data(i).Spectra,2);    
    labels((last+1):(last+nt)) = i;
    last = last+nt;
end

[~,n_pix] = size(hsi_data);

% classify by majority of K nearest neighbors
knn_out = zeros(n_pix,1);
idx = knnsearch(train',hsi_data','K',K);

for i=1:n_pix
   
    counts = zeros(n_class,1);
    for j=1:K
        counts(labels(idx(i,j))) = counts(labels(idx(i,j)))+1;
    end
    [~,max_i] = max(counts);
    knn_out(i) = max_i;
        
end

end
