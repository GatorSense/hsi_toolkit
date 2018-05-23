function [scores,det_out,names] = segsigdet_gauntlet(hsi,img,tgt_sig,segments,filt,title_str,roc_figh)
%
%function [scores,det_out,names] = segsigdet_gauntlet(hsi,img,tgt_sig,segments,filt,title_str,roc_figh)
%
% run the gauntlet of segmented signature detector tests for
% the given image and signature and segments
%
% 9/2012 - Taylor C. Glenn - tcg@cise.ufl.edu

if ~exist('roc_figh','var'); roc_figh = []; end

%-------------------------------------------------------------------------------
% Spectral angle mapper
i = 1;
det_out{i} = segmented(@sam_detector,img,tgt_sig,segments);
names{i} = 'Segmented Spectral Angle Mapper';

% Spectral matched filter
i = i+1;
det_out{i} = segmented(@smf_detector,img,tgt_sig,segments);
names{i} = 'Segmented Spectral Matched Filter';

% Constrained energy minimization
i = i+1;
det_out{i} = segmented(@cem_detector,img,tgt_sig,segments);
names{i} = 'Segmented Constrained Energy Minimization Filter';

% Adaptive coherence/cosine estimator (ACE)
i = i+1;
det_out{i} = segmented(@ace_detector,img,tgt_sig,segments);
names{i} = 'Segmented Adaptive coherence/cosine estimator (ACE)';

% Subpixel Spectral matched filter
i = i+1;
det_out{i} = segmented(@spsmf_detector,img,tgt_sig,segments);
names{i} = 'Segmented Subpixel Spectral Matched Filter';

% class conditional matched filter
i = i+1;
det_out{i} = segmented(@ccmf_detector,img,tgt_sig,segments,4);
names{i} = 'Segmented Class Conditional Matched Filter';

% Pairwise adaptive linear matched filter
i = i+1;
det_out{i} = segmented(@palm_detector,img,tgt_sig,segments,4);
names{i} = 'Segmented Pairwise Adaptive Linear Matched Filter';

% Orthogonal subspace projection
i = i+1;
det_out{i} = segmented(@osp_detector,img,tgt_sig,segments,2);
names{i} = 'Segmented Orthogonal Subspace Projection Filter';

n_det = i;

scores = cell(1,n_det);
for i=1:n_det
    scores{i} = score_hylid_perpixel(hsi,det_out{i},filt,names{i});
end

if ~isempty(roc_figh)
        
    figure(roc_figh);    
    PlotBullwinkleRoc(scores,title_str,names);
        
end
