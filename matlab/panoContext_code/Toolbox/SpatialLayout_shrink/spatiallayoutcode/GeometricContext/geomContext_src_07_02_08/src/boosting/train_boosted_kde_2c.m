function density = train_boosted_kde_2c(data, labels, ranges, num_iter)
% Try to learn ln(P(x1, x2 | +)/P(x1, x2 | -), where + indicates that a pair of points,
% x, are in the same cluster and - indicates that the pair are in different
% clusters.  Use boosting to estimate the parameters of the density in a
% naive structure.
% density(num_features).{x, log_ratio}


ndata = numel(labels);
nfeatures = size(data, 2);

if ~exist('ranges') || isempty(ranges)
    ranges = cell(nfeatures, 1);
end
    
for f = 1:numel(ranges)
    if isempty(ranges{f})
        ranges{f} = 'unbounded';
    end
end

pos_indices = find(labels==1);
neg_indices = find(labels==-1);
npos = length(pos_indices);
nneg = length(neg_indices);
%weights = 1/num_data*ones(num_data, 1);
% weights = zeros(ndata, 1);
% weights(pos_indices) = 1/2/npos;
% weights(neg_indices) = 1/2/nneg;
weights = repmat(1/ndata, [ndata 1]);
pos_data = data(pos_indices, :);
neg_data = data(neg_indices, :);

kernel_width = zeros(nfeatures, 1);

% get optimal kernel width and xrange
% this window parameter is recommended by Silverman(1986) for non-normal plots
for f = 1:nfeatures
         
    y = data(:, f);
    s = std(y);
    sorty = sort(y);
    n = length(y);
    iq = sorty(round(3/4*n))-sorty(round(n/4));
    density(f).kernelwidth = 0.9*min(s, iq/1.34)*(1/n)^(1/5); 
    %disp(num2str(density(f).kernelwidth))
    [tmp, density(f).x] = ksdensity(y, 'weights', weights, ...
        'width', density(f).kernelwidth, 'support', ranges{f});
    density(f).log_ratio = zeros(size(density(f).x));
end

total_confidences = zeros(ndata, 1);
for m = 1:num_iter
    pos_weights = weights(pos_indices);
    neg_weights = weights(neg_indices);
    tmp_confidences = zeros(ndata, 1);
    % update densities
    for f = 1:nfeatures       
        pos_f = ksdensity(pos_data(:, f), density(f).x, 'weights', pos_weights, ...
            'width', density(f).kernelwidth, 'support', ranges{f});
        pos_f = pos_f + 1/(ndata/2);
        pos_f = pos_f / sum(pos_f) * sum(pos_weights);        
        neg_f = ksdensity(neg_data(:, f), density(f).x, 'weights', neg_weights, ...
            'width', density(f).kernelwidth, 'support', ranges{f});       
        neg_f = neg_f + 1/(ndata/2);
        neg_f = neg_f / sum(neg_f) * sum(neg_weights);        
                                
        tmp_ratio = (log(pos_f)-log(neg_f));
        %density(f).log_ratio = density(f).log_ratio + tmp_ratio;
        curr_ratio(f).log_ratio = tmp_ratio;
        tmp_confidences = tmp_confidences + ...
            test_boosted_kde_2c(tmp_ratio, density(f).x, data(:, f))';
    end
    
    % get alpha parameter for confidence
    alpha = fminbnd(@compute_expected_confidence, 0.001, 2.0, [], labels, tmp_confidences, weights); 
    %disp(['alpha = ' num2str(alpha)]);
    for f = 1:nfeatures
        density(f).log_ratio = density(f).log_ratio + alpha*curr_ratio(f).log_ratio;
    end
    
    weights = weights .* exp(-alpha*tmp_confidences.*labels);
    sumw = sum(weights);
    %disp(['sum w = ' num2str(sumw)]);
    weights = weights / sumw;
    
    total_confidences = total_confidences + alpha*tmp_confidences;
    
    disp(['training error:  n_err = ' num2str(mean(total_confidences(neg_indices)>=0)) ' p_err = ' ...
            num2str(mean(total_confidences(pos_indices)<0))]);    
    if 0    
    [f1, x] = ksdensity(total_confidences(pos_indices));
    f2 = ksdensity(total_confidences(neg_indices), x);
    figure(1), plot(x, f1, 'b', x, f2, 'r');    
    pause(0.5);
    end

end




function expected_confidence = compute_expected_confidence(alpha, labels, confidences, weights)        
% used to choose alpha
% alpha is chosen so that the expected confidence:
% E(weight*exp(-alpha*label*confidence)*label*confidence) = 0
% Note: Absolute value used so that minimization results in best confidence
% confidence.  Normalization is unnecesary for finding alpha, though.
new_weights = exp(-alpha*labels.*confidences).*weights;    
expected_confidence = sum(new_weights);
%expected_confidence = abs(sum(new_weights.*labels.*confidences)/sum(new_weights));
    