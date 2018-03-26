function density = boostKernelDensity(X, y, numIter, w, posprior)
% Boosted naive kernel density classifier/regressor for two classes. 
% Attempts to estimate 1./(1+exp(-f(X))) = y where f(X) is the log ratio
% of class likelihoods.  Uses adaboost with confidence-weighted predictions
%
% INPUT
% X(ndata, nfeatures) - the features
% y(ndata, 1) - class confidence ranges from 0 to 1
% numIter - number of boosting iterations
% w(ndata, 1) - set of initial weights (optional)
% posprior - prior for positive class (optional)
%
% OUTPUT
% classifier.{x, log_ratio}


[ndata, nfeatures]  = size(X);

if ~exist('w')
    w = ones(ndata, 1);
end
if ~exist('posprior')
    posprior = mean(y);
end

% regression
regression = 0;
if any(y~=0 & y~=1)
    regression  = 1;
    X = repmat(X, 2, 1);
    w = [y.*w ; (1-y).*w];
    pind = [1:ndata];
    nind =[(ndata+1):(2*ndata)];
    ndata = ndata*2;
else
    pind = find(y==1);
    nind = find(y==0);    
end

w(pind) = w(pind) / sum(w(pind)) * posprior;
w(nind) = w(nind) / sum(w(nind)) * (1-posprior);

labels = zeros(ndata, 1);
labels(pind) = 1;
labels(nind) = -1;

kernel_width = zeros(nfeatures, 1);

% get optimal kernel width and xrange as recommended by Silverman(1986) 
% for non-normal data
for f = 1:nfeatures
    
    if regression
        n = ndata/2;
    else
        n = ndata;
    end
    s = std(X(:, f));
    sorty = sort(X(:, f));
    iq = sorty(round(3/4*n))-sorty(round(n/4));
    density(f).kernelwidth = 0.9*min(s, iq/1.34)*(1/n)^(1/5) * 4; % multiplying by 4 as hack
    [tmp, density(f).x] = ksdensityw(X(:, f), w, 'width', density(f).kernelwidth);
    [tmp, firstx] = min(density(f).x + (density(f).x<0)*1e10); % smallest x above 0
    density(f).x = density(f).x(firstx:end);        
    
    density(f).log_ratio = zeros(size(density(f).x))';
end

total_confidences = zeros(ndata, 1);
for m = 1:numIter
    tmp_confidences = zeros(ndata, 1);
    % update densities
    for f = 1:nfeatures
        pprob = ksdensityw(X(pind, f), w(pind), density(f).x, 'width', density(f).kernelwidth)';
        pprob = pprob + 1/ndata;
        pprob = pprob / sum(pprob) * sum(w(pind));
        nprob = ksdensityw(X(nind, f), w(nind), density(f).x, 'width', density(f).kernelwidth)';
        nprob = nprob + 1/ndata;
        nprob = nprob / sum(nprob) * sum(w(nind));                                        
        tmp_ratio = (log(pprob)-log(nprob));
        curr_ratio(f).log_ratio = tmp_ratio;
        tmp_confidences = tmp_confidences + get_likelihood_ks(tmp_ratio, density(f).x, X(:, f));
    end
    
    % get alpha parameter for confidence
    alpha = fminbnd(@compute_expected_confidence, 0.01, 2.0, [], labels, tmp_confidences, w); 
    for f = 1:nfeatures
        density(f).log_ratio = density(f).log_ratio + alpha*curr_ratio(f).log_ratio;
    end
    disp(num2str(alpha))
    total_confidences = total_confidences + alpha*tmp_confidences;    
    
    %w = 1 ./ (1+exp(labels.*total_confidences));        
    %w = w / sum(w);      
    
    w = w .* exp(-alpha*tmp_confidences.*labels);
    sumw = sum(w);
    disp(num2str(sumw))
    
    w = w / sumw;
    
    if regression
        disp(['confidence error: ' num2str(mean(...
            abs(1./(1+exp(-total_confidences(1:length(y))))-y)))]);
    else
        disp(['training error:  n_err = ' num2str(mean(total_confidences(nind)>=0)) ' p_err = ' ...
                num2str(mean(total_confidences(pind)<0)) 'total: ' ...
                num2str(mean(total_confidences.*labels < 0))]);    
    end
    if 0    
    [f1, x] = ksdensity(total_confidences(pos_indices));
    f2 = ksdensity(total_confidences(neg_indices), x);
    figure(1), plot(x, f1, 'b', x, f2, 'r');    
    pause(0.5);
    end
    pause(0.1);
end
    
