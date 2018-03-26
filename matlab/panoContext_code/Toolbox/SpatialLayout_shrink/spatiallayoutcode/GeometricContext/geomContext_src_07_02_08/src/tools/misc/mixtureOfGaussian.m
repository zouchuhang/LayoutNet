function gmm = mixtureOfGaussian(data, K, niter, fullCov, varargin)
% Mixture of K gaussians
% data(ndata, nvariables) - the data
% K - the number of gaussians in the mixture
% niter - the number of iterations
% fullCov - 1 if a full covariance matrix should be used, otherwise
%                   a diagonal matrix will be used
% varargin{1} - weights for the data (optional)
% gmm.mu{K} - the means of the mixtures
% gmm.sigma{K} - the variance matrices
% gmm.priors(K) - the mixture priors

[ndata, ndim] = size(data);

if ndata == 0
    disp('warning no data')
    gmm.mu = {};
    gmm.sigma = {};
    gmm.priors = [];
    return;
end
   

if length(varargin)==0
    w = ones(ndata, 1);
else
    w = varargin{1}(:);
end

w = w / sum(w) * ndata;

cumw = cumsum(w / sum(w));

% initialize mu by randomly selected data values 
r = rand(K, 1);
for k = 1:K
    for i = 1:ndata
        if cumw(i)>=r(k)
            mu{k} = data(i, :);
            break;
        end
    end
end

% initialize all covariances to be covariance of data
%if fullCov
%    sigma(1:K) = {cov(data)};
%else
%    sigma(1:K) = {diag(var(data))};
%end
sigma(1:K) = {diag(ones(1, ndim))};

priors = ones(1, K) / K;

for iter = 1:niter
        
    
    % get posteriors: P(k | data) = P(data|k)P(k) / sum_k(P(data|k))
    pC = zeros(ndata, K);   
    for k = 1:K
        invSigma = inv(sigma{k});
        for i = 1:ndata
            pC(i, k) = 1 / (sqrt(2*pi)*sqrt(det(sigma{k}))) * ...
                exp(-1/2 * (data(i, :)-mu{k}) * invSigma * (data(i, :)-mu{k})') ...
                * priors(k);
        end
    end
    
    if 0 
    for k = 1:K
        disp(['k: ' num2str(k)])
        disp(['mu: ' num2str(mu{k})])
        disp(['sigma: ' num2str(det(sigma{k}))])
        disp(['prior: ' num2str(priors(k))])        
        disp(['min: ' num2str(min(pC(:, k)))])
        disp(['max: ' num2str(max(pC(:, k)))])
    end    
    end

    %disp(num2str(priors))
    
    pC = pC + eps; % avoid divide by zero error
    
    sumPc = sum(pC, 2);
    for k = 1:K
        pC(:, k) = pC(:, k) ./ sumPc .* w;
    end    
    kweight = sum(pC, 1);    
    
    priors = kweight / sum(kweight);
    
    % mean estimation    
    for k = 1:K
        sumkX = pC(:, k)' * data;
        mu{k} = sumkX / kweight(k);
    end
    
    % covariance matrix estimation
    if fullCov
        for k = 1:K        
            diffs = data - (ones(ndata, 1) * mu{k});
            diffs = diffs .* (sqrt(pC(:,k))*ones(1, ndim));
            sigma{k} = (diffs' * diffs) / kweight(k);
        end
    else
        for k = 1:K
            diffs = data - (ones(ndata, 1) * mu{k});
            sigma{k} = diag(sum((diffs.*diffs).*(pC(:, k)*ones(1, ndim)), 1) ./ kweight(k));
        end
    end
end
  
gmm.mu = mu;
gmm.sigma = sigma;
gmm.priors = priors;

