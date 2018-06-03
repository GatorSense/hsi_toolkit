function [mergedData] = dimReduction(img, Parameters)

%function [dimReductionStruct] = dimReduction(img, Parameters)
%
% Hierarchical Dimensionality Reduction
% Computes KL-divergence between every pair of bands in image and merges
% bands that are similar based on KL-divergence.
%
% Inputs:  
%   img: hyperspectral data cube (n_row x n_col x n_bands)
%   Parameters: parameter structure defined by dimReductionParameters.m
%
% Author: Alina Zare
% Email Address: azare@ufl.edu
% Created: September 12, 2008
% Latest Revision: June 3, 2018
%%

[numRows, numCols, numDims] = size(img);

if(nargin < 2)
    Parameters = dimReductionParameters();
end

%Perform Hierarchical Dimensionality Reduction
type = Parameters.type; %Type of Hierarchy
showH = Parameters.showH; %Set to 1 to show clustering, 0 otherwise
maxNumClusters = Parameters.numBands; 
NumCenters = Parameters.NumCenters; %Number of centers used in computing KL-divergence

InputData = reshape(img, [numRows*numCols, numDims]);
[~, KLDivergencesList, ~] = computeKLDivergencesBetweenBands(InputData, NumCenters);

Hierarchy = linkage(KLDivergencesList, type);
if(showH)
    D = dendrogram(Hierarchy, 0);
end
band_clusters = cluster(Hierarchy,'maxclust',maxNumClusters) ;

mergedData = zeros(maxNumClusters, numRows*numCols);
for i = 1:maxNumClusters
    mergedData(i,:) = mean(InputData(:,band_clusters == i),2);
end

mergedData = reshape(mergedData', [numRows, numCols, maxNumClusters]);



end

%%
function [KLDivergences, KLDivergencesList, hists] = computeKLDivergencesBetweenBands(InputData, NumCenters)


DataList = InputData/max(max(InputData));

%Compute Histograms
Centers = [1/NumCenters:1/NumCenters:1];
hists = hist(DataList', Centers);
hists = hists+eps;

%Compute KL-divergences
KLDivergences = zeros([size(InputData,2), size(InputData,2)]);
for i = 1:size(DataList,2)
    for j = 1:size(DataList,2)
        KLDivergences(i,j) = sum(hists(i,:).*log(hists(i,:)./(hists(j,:)))) + sum(hists(j,:).*log(hists(j,:)./(hists(i,:))));
    end
end

%Sort in List order for linkage algorithm
temp = KLDivergences - diag(diag(KLDivergences));
KLDivergencesList = squareform(temp);

end

