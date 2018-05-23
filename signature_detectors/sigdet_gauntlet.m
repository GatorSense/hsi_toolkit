function [scores,det_out,names] = sigdet_gauntlet(hsi,img,tgt_sig,filt,title_str,roc_figh)
%
%function [scores,det_out,names] = sigdet_gauntlet(hsi,img,tgt_sig,filt,title_str,roc_figh)
%
% run the gauntlet of signature detector tests for
% the given image and signature
%

if ~exist('roc_figh','var'); roc_figh = []; end

%-------------------------------------------------------------------------------
% Spectral angle mapper
i = 1;
det_out{i} = sam_detector(img,tgt_sig);
names{i} = 'Spectral Angle Mapper';

% Spectral matched filter
i = i+1;
det_out{i} = smf_detector(img,tgt_sig);
names{i} = 'Spectral Matched Filter';

% Constrained energy minimization
i = i+1;
det_out{i} = cem_detector(img,tgt_sig);
names{i} = 'Constrained Energy Minimization Filter';

% Adaptive coherence/cosine estimator (ACE)
i = i+1;
det_out{i} = ace_detector(img,tgt_sig);
names{i} = 'Adaptive coherence/cosine estimator (ACE)';

% Subpixel Spectral matched filter
i = i+1;
det_out{i} = spsmf_detector(img,tgt_sig);
names{i} = 'Subpixel Spectral Matched Filter';

% class conditional matched filter
i = i+1;
det_out{i} = ccmf_detector(img,tgt_sig,[],8);
names{i} = 'Class Conditional Matched Filter';

% Pairwise adaptive linear matched filter
i = i+1;
det_out{i} = palm_detector(img,tgt_sig,[],8);
names{i} = 'Pairwise Adaptive Linear Matched Filter';

% Orthogonal subspace projection
i = i+1;
det_out{i} = osp_detector(img,tgt_sig,[],10);
names{i} = 'Orthogonal Subspace Projection Filter';

% %-------------------------------------------------------------------------------
% % Detectors and False Alarm Mitigation from Eismann 14.4
% 
% % Quadratic Matched Filter
% det_out{i} = qmf_detector(img,tgt_sig,0.1*eye(n_dim_pca));
% 
% score_hylid_perpixel(hsi,qmf_out,0,filt,'Quadratic Matched Filter',figure(26),figure(27));
% 
% 
% % False Alarm Mitigation by Subpixel Replacement Model
% fam_out = fam_statistic(img,tgt_sig);
% smf_fam = smf_out;
% smf_fam(fam_out > 100) = 0;
% 
% score_hylid_perpixel(hsi,smf_fam,filt,'Spectral Matched Filter w/ False Alarm Mitigation',figure(28),figure(29));
% 
% % Mixture Tuned Matched Filter (infeasibility metric)
% % fixme: needs work
% [mtmf_out,mtmf_alpha] = mtmf_statistic(hsi_img,target_sig);
% smf_mtmf = mtmf_alpha*1e4;
% smf_mtmf(log10(mtmf_out+1) > 8) = 0;
% 
% score_hylid_perpixel(hsi,smf_mtmf,filt,'Mixture Tuned Matched Filter',figure(30),figure(31));
% 
% % Finite Target Matched Filter
% % fixme: needs work
% ftmf_out = ftmf_detector(img,tgt_sig,0.25);
% 
% score_hylid_perpixel(hsi,ftmf_out,filt,'Finite Target Matched Filter',figure(32),figure(33));
% 
% % Least angle regression
% [lar_out] = lar_statistic(hsi_img,target_sig,smf_out > 1.5);
% smf_lar = smf_out;
% smf_lar(lar_out < 1.5) = 0;
% 
% score_hylid_perpixel(hsi,smf_lar,filt,'Least Angle Regression with Matched Filter',figure(34),figure(35));

n_det = i;

scores = cell(1,n_det);
for i=1:n_det
    scores{i} = score_hylid_perpixel(hsi,det_out{i},filt,names{i});
end

if ~isempty(roc_figh)
        
    figure(roc_figh);    
    PlotBullwinkleRoc(scores,title_str,names);
        
end


