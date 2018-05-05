function det_out = anomaly_det_demo(hsi_img, mask, wavelengths)
%
%det_out = anomaly_det_demo(hsi_img, mask)
%
% Demo that runs all anomaly detectors in the hsi_toolkit
% 
% How to run using example sub MUUFL Gulfport data:
%	(1) load an_hsi_img_sub_for_demo.mat file 
%	(2) run anomaly_det_demo(hsi_img_sub, mask_sub, wavelengths)
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%  wavelengths - 1 x n_band vector listing wavelength values for hsi_img in nm
%
% outputs:
%   det_out - cell array of detector output images
%
% 5/5/2018 - Alina Zare

det_out = {};

[beta_out] = beta_anomaly(hsi_img,mask);
det_out{end+1}.result = beta_out; 
det_out{end}.method = 'beta\_anomaly';

det_out{end+1}.parameters.n_cluster = 8;
[cbad_out,cluster_img] = cbad_anomaly(hsi_img,mask,det_out{end}.parameters.n_cluster);
det_out{end}.result = cbad_out; 
det_out{end}.cluster_img = cluster_img;
det_out{end}.method = 'cbad\_anomaly';

det_out{end+1}.parameters.n_dim_bg = 3;
det_out{end}.parameters.n_dim_tgt = [];
det_out{end}.parameters.tgt_orth = true;
[csd_out] = csd_anomaly(hsi_img,det_out{end}.parameters.n_dim_bg,det_out{end}.parameters.n_dim_tgt,det_out{end}.parameters.tgt_orth);
det_out{end}.result = csd_out; 
det_out{end}.method = 'csd\_anomaly';

det_out{end+1}.parameters.n_cluster = 8;
[fcbad_out,cluster_img] = fcbad_anomaly(hsi_img,mask,det_out{end}.parameters.n_cluster);
det_out{end}.result = fcbad_out; 
det_out{end}.cluster_img = cluster_img;
det_out{end}.method = 'fcbad\_anomaly';

det_out{end+1}.parameters.n_comp = 8;
[gmm_out] = gmm_anomaly(hsi_img,mask,det_out{end}.parameters.n_comp);
det_out{end}.result = gmm_out; 
det_out{end}.method = 'gmm\_anomaly';

det_out{end+1}.parameters.n_comp = 8;
[gmrx_out] = gmrx_anomaly(hsi_img,mask,det_out{end}.parameters.n_comp );
det_out{end}.result = gmrx_out; 
det_out{end}.method = 'gmrx\_anomaly';

[md_out] = md_anomaly(hsi_img,mask);
det_out{end+1}.result = md_out; 
det_out{end}.method = 'md\_anomaly';

det_out{end+1}.parameters.guard_win = 3;
det_out{end}.parameters.bg_win = 3;
[rx_out] = rx_anomaly(hsi_img,mask,det_out{end}.parameters.guard_win,det_out{end}.parameters.bg_win);
det_out{end}.result = rx_out; 
det_out{end}.method = 'rx\_anomaly';

det_out{end+1}.parameters.n_dim_ss = 3;
det_out{end}.parameters.guard_win = 3;
det_out{end}.parameters.bg_win = 3;
[ssrx_out] = ssrx_anomaly(hsi_img, det_out{end}.parameters.n_dim_ss, det_out{end}.parameters.guard_win, det_out{end}.parameters.bg_win);
det_out{end}.result = ssrx_out; 
det_out{end}.method = 'ssrx\_anomaly';

figure;
numR = 4;
numC = 3; 
subplot(numR,numC,1); imagesc(get_RGB(hsi_img, wavelengths)); title('RGB image');
subplot(numR,numC,2); imagesc(mask); title('valid mask');
for i = 1:length(det_out)
	subplot(numR, numC, i+2); imagesc(det_out{i}.result); title(det_out{i}.method);
end