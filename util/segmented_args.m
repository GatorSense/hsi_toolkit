function det_out = segmented_args(detector_fn,hsi_img,tgt_sig,segments,varargin)
%
%function det_out = segmented_args(detector_fn,hsi_img,tgt_sig,segments,varargin)
%
% Segmented Detector Wrapper
%  uses any detector with the signature detector(img,tgt_sig,mask,args...)
%  as a segmented detector over the given segments
%  Distributes cell array arguments from varargin to corresponding segment
%
% inputs:
%  detector_fn - function handle for wrapped detector
%  hsi_img - n_row x n_col x n_band hyperspectral image
%  tgt_sig - target signature (n_band x 1 - column vector)
%  segments - cell array of segment masks, n_row x n_col binary images
%  varargin - variable array of arguments passed to the detector function
%             cell array arguments are distributed to corresponding segment
%
% outputs:
%  det_out - detector output image, concatenation of outputs from each segment
%            NaN valued in pixels not contained by a segment
%
% 8/20/2012 - Taylor C. Glenn 
%

[n_row,n_col,~] = size(hsi_img);

n_seg = numel(segments);

det_out = NaN(n_row,n_col);

for i=1:n_seg
    
    vararg = cell(1,numel(varargin));
    for j=1:numel(varargin)
        if iscell(varargin{j})
            vararg{j} = varargin{j}{i};
        else
            vararg{j} = varargin{j};
        end
    end
    
    seg_out = detector_fn(hsi_img,tgt_sig,segments{i},vararg{:});
    
    det_out(segments{i}) = seg_out(segments{i});
    
end

end