function [det_out,varargout] = segmented(detector_fn,hsi_img,tgt_sig,segments,varargin)
%
%function [det_out,varargout] = segmented(detector_fn,hsi_img,tgt_sig,segments,varargin)
%
% Segmented Detector Wrapper
%  uses any detector with the signature detector(img,tgt_sig,mask,args...)
%  as a segmented detector over the given segments
%
% inputs:
%  detector_fn - function handle for wrapped detector
%  hsi_img - n_row x n_col x n_band hyperspectral image
%  tgt_sig - target signature (n_band x 1 - column vector)
%  segments - cell array of segment masks, n_row x n_col binary images
%  varargin - variable array of arguments passed to the detector function
%
% outputs:
%  det_out - detector output image, concatenation of outputs from each segment
%            NaN valued in pixels not contained by a segment
%  varargout - other outputs, assumed to be images
%
% 8/20/2012 - Taylor C. Glenn 
%

[n_row,n_col,~] = size(hsi_img);

n_seg = numel(segments);

det_out = NaN(n_row,n_col);        

for i=1:n_seg

    if(isempty(varargin))
        seg_out = detector_fn(hsi_img,tgt_sig,segments{i});
    else
        seg_out = detector_fn(hsi_img,tgt_sig,segments{i},varargin{:});
    end

    det_out(logical(segments{i})) = seg_out(logical(segments{i}));

end


end