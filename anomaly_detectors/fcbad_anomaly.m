function [fcbad_out,cluster_img] = fcbad_anomaly(hsi_img,mask,n_cluster)
%
%fcbad_out = fcbad_anomaly(hsi_img,mask,n_cluster)
%
% Fuzzy Cluster Based Anomaly Detection (FCBAD)
% Ref: Hytla, Patrick C., et al. "Anomaly detection in hyperspectral imagery: comparison of methods using diurnal and seasonal data." Journal of Applied Remote Sensing 3.1 (2009): 033546
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%  n_cluster - number of clusters to use
%
% outputs:
%  fcbad_out - detector output image
%  cluster_img - cluster label image
%
% 8/8/2012 - Taylor C. Glenn 
% 5/5/2018 - Edited by Alina Zare
    
    addpath(fullfile('..','util'));
    [fcbad_out,cluster_img] = img_det(@fcbad_det,hsi_img,[],mask,n_cluster);
end

function [fcbad_data,idx] = fcbad_det(hsi_data,~,n_cluster)
    [n_band,n_pix] = size(hsi_data);

    % cluster the data

    %options(1): exponent for the partition matrix U (default: 2.0)
    %options(2): maximum number of iterations (default: 100)
    %options(3): minimum amount of improvement (default: 1e-5)
    %options(4): info display during iteration (default: 1)

    opt = [2 500 1e-6 0];

    [C,U] = fcm(hsi_data',n_cluster,opt);  %this requires the fuzzy logic toolbox

    % make a crisp partitioning
    idx = max(U,[],1);

    % get cluster statistics
    mu = cell(1,n_cluster);
    siginv = cell(1,n_cluster);
    for i=1:n_cluster
        
        mu{i} = C(i,:)';

        % computer a membership weighted covariance for the cluster
        sigma = zeros(n_band,n_band);
        for j=1:n_pix
           
            z = hsi_data(:,j) - mu{i};
            
            sigma = sigma + U(i,j) * z*z';    
        end
        sigma = sigma/sum(U(i,:));
        
        siginv{i} = pinv(sigma);
    end

    % compute total membership weighted Mahalanobis distances
    fcbad_data = zeros(1,n_pix);

    for j=1:n_pix
        m_dists = zeros(n_cluster,1);
        for i=1:n_cluster
            
            z = hsi_data(:,j) - mu{i};
            m_dists(i) = z'*siginv{i}*z;
        end
        
        fcbad_data(j) = U(:,j)'*m_dists;
    end
end