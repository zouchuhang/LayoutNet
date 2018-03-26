function px = gmmEval(x, gmm)
% Computes the probability of the data given a Guassian mixture model
% x(ndata, ndim) - the data to be evaluated
% gmm.mu{K}, gmm.sigma{K}, gmm.priors{K} - the mixture parameters

K = length(gmm.priors);
[ndata, ndim] = size(x);

px = zeros(ndata, 1);

if K ==0
    return;
end




for k = 1:K
    invSigma = inv(gmm.sigma{k});
    for i = 1:ndata
        px(i) = px(i) + 1/(sqrt(2*pi)*sqrt(det(gmm.sigma{k}))) * ...
            exp(-1/2 * (x(i, :)-gmm.mu{k}) * invSigma * (x(i, :)-gmm.mu{k})') ...
            * gmm.priors(k);
    end
end

