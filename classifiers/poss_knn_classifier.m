function pknn_img = poss_knn_classifier(hsi_img,train_data,K,m,eta)
%function pknn_img = poss_knn_classifier(hsi_img,train_data,mask,k,m,eta)
%
% Possibilistic K nearest nieghbors classifier
%
% Ref: Frigui, Hichem, and Paul Gader. "Detection and discrimination of land mines in ground-penetrating radar based on edge histogram descriptors and a possibilistic $ k $-nearest neighbor classifier." IEEE Transactions on Fuzzy Systems 17.1 (2009): 185-199.
%
% Inputs:
%   hsi_img: hyperspectral data cube (n_rows x n_cols x n_bands)
%   train_data: structure containing training data 
%       train_data(i).Spectra: matrix containing training data from class i
%   mask: binary image indicating where to apply classifier
%   K:  number of neighbors to use during classification
%   m:  fuzzifier (usually = 2)
%   eta: eta parameter to determine what is an outlier
%
% Outputs:
%   pknn_img: class membership matrix (n_row x n_col x n_class)
%
% 6/3/2018 - Alina Zare

[n_rows, n_cols, n_band] = size(hsi_img);
hsi_data = reshape(hsi_img, [n_rows*n_cols, n_band])';

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

%compute mu weights for training data
[idx_train] = knnsearch(train',train','K',K);
idx_labels = labels(idx_train);
mu = zeros(n_train, n_class);
for i=1:n_train
   u_lab = unique(idx_labels(i,:));
   for j = 1:length(u_lab)
       if(u_lab(j) == labels(i))
          count = sum(idx_labels(i,:) == u_lab(j));
          mu(i,u_lab(j)) = 0.51 + (count/K)*0.49;
       else
          count = sum(idx_labels(i,:) == u_lab(j));
          mu(i,u_lab(j)) = (count/K)*0.49;
       end
   end
end
  
% classify by using weighted K nearest neighbors
pknn_out = zeros(n_pix,n_class);
[idx, D] = knnsearch(train',hsi_data','K',K);
weights = D-eta;
weights(weights < 0) = 0;
weights = 1./(1+(weights.^(2/(m-1))+eps));
for i=1:n_pix
    pknn_out(i,:) = sum(repmat(weights(i,:)', [1, n_class]).*mu(idx(i,:),:)); 
end

pknn_img = zeros(n_rows, n_cols, n_class);
for i = 1:n_class
    pknn_img(:,:,i) = reshape(pknn_out(:,i), [n_rows, n_cols]);
end

end
