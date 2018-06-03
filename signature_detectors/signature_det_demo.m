function det_out = signature_det_demo(img, tgt_sig, mask, wavelengths, gtImage)
%
% det_out = signature_det_demo(img, tgt_sig, mask, wavelengths)
% 
% Demo that runs all signature detectors in the hsi_toolkit
% 
% How to run using example sub MUUFL Gulfport data:
% 	(1) load an_hsi_img_for_tgt_det_demo.mat file 
% 	(2) run signature_det_demo(hsi_sub, tgt_spectra, [], wavelengths, gtImg_sub)
% 
% inputs:
%  img - n_row x n_col x n_band hyperspectral image
%  tgt_sig - n_band x 1 target signature vector
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%  wavelengths - 1 x n_band vector listing wavelength values for hsi_img in nm
% 
% outputs:
%   det_out - cell array of detector outputs
% 
% 6/2/2018 - Alina Zare

addpath(fullfile('..','util'));

%% Signature Detectors

% Spectral angle mapper
i = 1;
det_out{i}.result = sam_detector(img,tgt_sig);
det_out{i}.method = 'Spectral Angle Mapper';

% Spectral matched filter
i = i+1;
det_out{i}.result = smf_detector(img,tgt_sig);
det_out{i}.method = 'Spectral Matched Filter';

% Constrained energy minimization
i = i+1;
det_out{i}.result = cem_detector(img,tgt_sig);
det_out{i}.method = 'Constrained Energy Minimization Filter';
 
% Adaptive coherence/cosine estimator (ACE)
i = i+1;
det_out{i}.result = ace_detector(img,tgt_sig);
det_out{i}.method = 'Adaptive coherence/cosine estimator squared';
 
% Class conditional matched filter
i = i+1;
det_out{i}.params.nComp = 2;
det_out{i}.result = ccmf_detector(img,tgt_sig,mask,det_out{i}.params.nComp);
det_out{i}.method = 'Class Conditional Matched Filter';

% Pairwise adaptive linear matched filter
i = i+1;
det_out{i}.params.nComp = 2;
det_out{i}.result = palm_detector(img,tgt_sig,mask,det_out{i}.params.nComp);
det_out{i}.method = 'Pairwise Adaptive Linear Matched Filter';

% Orthogonal subspace projection
i = i+1;
det_out{i}.params.n_dim_ss = 10;
det_out{i}.result = osp_detector(img,tgt_sig,mask,det_out{i}.params.n_dim_ss);
det_out{i}.method = 'Orthogonal Subspace Projection Filter';

% Target Abundance Map
i = i+1;
det_out{i}.params.ems = squeeze(img(1:3,2,:))'; %need to provide background endmembers (can get them using SPICE unmixing)
det_out{i}.result = abd_detector(img,tgt_sig,mask,det_out{i}.params.ems);
det_out{i}.method = 'Target Abundance Map';

% ACE detector
i = i+1;
det_out{i}.result = acert_detector(img,tgt_sig,mask);
det_out{i}.method = 'Adaptive coherence/cosine estimator';

% ACE RX detector
i = i+1;  
det_out{i}.params.guard_win = 1;
det_out{i}.params.bg_win = 3;
det_out{i}.params.beta = 0.001;
det_out{i}.result = ace_rx_detector(img,tgt_sig,[],det_out{i}.params.guard_win,det_out{i}.params.bg_win,det_out{i}.params.beta);
det_out{i}.method = 'Adaptive coherence RX estimator';

% ACE SS detector
i = i+1;
det_out{i}.result = ace_ss_detector(img,tgt_sig,mask);
det_out{i}.method = 'Adaptive coherence/cosine estimator subspace';

% Adaptive Matched Subspace Detector
i = i+1;
det_out{i}.params.n_dim_tgt = 1;
det_out{i}.params.n_dim_bg = 3;
det_out{i}.result = amsd_detector(img,tgt_sig,mask,det_out{i}.params.n_dim_tgt,det_out{i}.params.n_dim_bg);
det_out{i}.method = 'Adaptive Matched Subspace Detector';

% Cluster Tuned Matched Filter
i = i+1;
det_out{i}.params.nCluster = 2;
det_out{i}.result = ctmf_detector(img,tgt_sig,det_out{i}.params.nCluster);
det_out{i}.method = 'Cluster Tuned Matched Filter';

% False Alarm Mitigation Statistic
i = i+1;
det_out{i}.result = fam_statistic(img,tgt_sig);
det_out{i}.method = 'False Alarm Mitigation Statistic';

% Hybrid Abundance Detector
i = i+1;
det_out{i}.params.ems = double(squeeze(img(1:3,2,:))'); %need to provide background endmembers (can get them using SPICE unmixing)
det_out{i}.params.n_comp = 2;
det_out{i}.result = ha_detector(img,double(tgt_sig),[],det_out{i}.params.ems,det_out{i}.params.n_comp);
det_out{i}.method = 'Hybrid Abundance Detector';

% Hybrid Structured Detector
i = i+1;
det_out{i}.params.ems = double(squeeze(img(1:3,2,:))'); %need to provide background endmembers (can get them using SPICE unmixing)
det_out{i}.result = hsd_detector(img,double(tgt_sig),[],det_out{i}.params.ems);
det_out{i}.method = 'Hybrid Structured Detector';

% Hybrid Structured RX Detector
i = i+1;
det_out{i}.params.ems = double(squeeze(img(1:3,2,:))'); %need to provide background endmembers (can get them using SPICE unmixing)
det_out{i}.params.guard_win = 1;
det_out{i}.params.bg_win = 3;
det_out{i}.params.beta = 0.001;
det_out{i}.result = hsd_rx_detector(double(img),double(tgt_sig),[],det_out{i}.params.guard_win,det_out{i}.params.bg_win,double(det_out{i}.params.ems),det_out{i}.params.beta);
det_out{i}.method = 'Hybrid Structured RX Detector';

% Hybrid Unstructured Detector
i = i+1;
det_out{i}.params.ems = double(squeeze(img(1:3,2,:))'); %need to provide background endmembers (can get them using SPICE unmixing)
det_out{i}.params.n_comp = 2;
det_out{i}.result = hua_detector(img,double(tgt_sig),[],det_out{i}.params.ems,det_out{i}.params.n_comp);
det_out{i}.method = 'Hybrid Unstructured Detector';

% Quadratic Spectral Matched Filter
i = i+1;
det_out{i}.params.tgt_cov = 0.1*eye(size(img,3));
det_out{i}.result = qmf_detector(img,tgt_sig,det_out{i}.params.tgt_cov);
det_out{i}.method = 'Quadratic Spectral Matched Filter';

% Spectral Matched Filter RX Detector
i = i+1;
det_out{i}.params.ems = double(squeeze(img(1:3,2,:))'); %need to provide background endmembers (can get them using SPICE unmixing)
det_out{i}.params.guard_win = 1;
det_out{i}.params.bg_win = 3;
det_out{i}.params.beta = 0.001;
det_out{i}.result = smf_rx_detector(img,tgt_sig,[],det_out{i}.params.guard_win,det_out{i}.params.bg_win);
det_out{i}.method = 'Spectral Matched Filter RX Detector';

%% Segmented Detector Examples

%Get Segments (using K-means here, but better ways to do this in general, see context-dependent methods for detection)
numK = 3;
[idx] = kmeans(reshape(img, [size(img,1)*size(img,2), size(img,3)]), numK);
idxIm = reshape(idx, [size(img,1), size(img,2)]);
for j = 1:numK
    segments{j} = zeros([size(img,1), size(img,2)]);
    segments{j}(idxIm == j) = 1;
end

% Segmented Spectral Angle Mapper
i = i+1;
det_out{i}.result = segmented(@sam_detector,img,tgt_sig,segments);
det_out{i}.method  = 'Segmented Spectral Angle Mapper';

% Segmented ACE
i = i+1;
det_out{i}.result = segmented(@ace_detector,img,tgt_sig,segments);
det_out{i}.method  = 'Segmented ACE';
    

%%  Visualize Results

figure;
numR = 4;
numC = ceil((i+1)/numR); 
subplot(numR,numC,1); imagesc(get_RGB(img, wavelengths)); title('RGB image');
subplot(numR,numC,2); imagesc(gtImage); title('Groundtruth image');
for i = 1:length(det_out)
	subplot(numR, numC, i+2); imagesc(det_out{i}.result); title(det_out{i}.method);
end

