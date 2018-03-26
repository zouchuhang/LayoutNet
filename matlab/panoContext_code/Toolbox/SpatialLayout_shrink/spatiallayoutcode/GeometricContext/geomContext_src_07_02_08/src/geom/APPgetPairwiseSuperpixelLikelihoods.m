function Z = APPgetPairwiseSuperpixelLikelihoods(features, density)
% Z = APPgetPairwiseSuperpixelLikelihoods(features, density)
% Gets the matrix of pairwise sp likelihoods

nElements = size(features, 1);
nFeatures = size(features, 2);
nDist = nElements*(nElements-1)/2;

% Y contains all pairwise likelihoods of being in same cluster
Y = zeros(nDist, 1);
count = 0;
for f = 1:nFeatures
    tmpY = pdist(features(:, f), 'cityblock');
    % minus used because log ratio is of P(x|+)/P(x|-), but we want inverse    
    %% NOW plus used for similarity
    Y = Y + getKdeLikelihood(density(f).log_ratio, density(f).x, tmpY);
end

Z = squareform(Y');
Z = Z.*(ones(nElements, nElements)-eye(nElements));