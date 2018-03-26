function classifier = train_boosted_dt_mc(features, cat_features, labels, ...
    num_iterations, num_nodes, stopval, init_weights, varargin)
%
%classifier = train_boosted_dt_mc(features, cat_features, labels, ...
%    num_iterations, num_nodes, stopval, init_weights, varargin)
%
% Train a classifier based on boosted decision trees.  Boosting done by the
% logistic regression version of Adaboost (Adaboost.L - Collins, Schapire,
% Singer 2002).  At each
% iteration, a set of decision trees is created for each class, with
% confidences equal to 1/2*ln(P+/P-) for that class, according to the
% weighted distribution.  Final classification is based on the largest
% confidence label (possibly incorporating a prior as h0(c) =
% 1/2*ln(Pc/(1-Pc)).  Weights are assigned as
% w(i,j) = 1 / (1+exp(sum{t in iterations}[yij*ht(xi, j)])).  

if length(varargin) == 1  % class names supplied
    gn = varargin{1};
    gid = zeros(size(labels));
    for c = 1:length(gn)
        ind = find(strcmp(labels, gn{c}));
        gid(ind) = c;
        if ~isempty(init_weights)
            disp([gn{c} ': ' num2str(sum(init_weights(ind)))]);
        else
            disp([gn{c} ': ' num2str(length(ind))]);
        end
    end
    ind = find(gid==0);
    gid(ind) = [];
    labels(ind) = [];
    features(ind, :) = [];
else    
    [gid, gn] = grp2idx(labels);    
end

if ~exist('stopval', 'var') || isempty(stopval)
    stopval = 0;
end
if ~exist('init_weights', 'var') 
    init_weights = [];
end

classifier.names = gn;

num_classes = length(gn);
num_data = length(gid);

if isempty(init_weights)
    init_weights = ones(num_data, 1)/num_data;
else
    init_weights = init_weights / sum(init_weights);
end

% if no examples from a class are present, create one dummy example for
% that class with very small weight
for c = 1:numel(gn)
    if ~any(gid==c)
        disp(['warning: no examples from class ' gn(c)])
        gid(end+1) = c;
        features(end+1, :) = zeros(size(features(end, 1)));
        num_data = num_data + 1;
        init_weights(end+1) = min(init_weights)/2;        
    end
end

all_conf = zeros(num_data, num_classes);
for c = 1:num_classes

    disp(['class: ' num2str(gn{c})]);    
    y = (gid == c)*2-1;
    cl = [-1 1];
    nc = 2;
    w = zeros(num_data, 1);
    cw = zeros(num_classes, 1);  
    for i = 1:2
        indices = find(y==cl(i));
        %count = sum(init_weights(indices));
        %w(indices) = init_weights(indices) / count / 2;
        w(indices) = init_weights(indices);
        
        if cl(i)==1
            %classifier.h0(c) = log(count / (1-count));
            classifier.h0(c) = 0;
        end
        
    end
        
    data_confidences = zeros(num_data, 1);
    aveconf = [];
    
    for t = 1:num_iterations
        % learn decision tree based on weighted distribution
        dt = treefitw(features, y, w, 1/num_data/2, 'catidx', cat_features, 'method', 'classification', 'maxnodes', num_nodes*4);
        [tmp, level] = min(abs(dt.ntermnodes-num_nodes));
        dt = treeprune(dt, 'level', level-1);

        % assign partition confidences
        pi = (strcmp(dt.classname{1},'1')) + (2*strcmp(dt.classname{2},'1'));
        ni = (strcmp(dt.classname{1},'-1')) + (2*strcmp(dt.classname{2},'-1'));
        
        classprob = dt.classprob;
        confidences = 1/2*(log(classprob(:, pi)) - log(classprob(:, ni)));             

        % assign weights
        [class_indices, nodes, classes] = treeval(dt, features);        
        data_confidences = data_confidences + confidences(nodes);
        
        w = 1 ./ (1+exp(y.*data_confidences));        
        w = w / sum(w);   
                
        %disp(['c: ' num2str(mean(1 ./ (1+exp(-y.*data_confidences)))) '  e: ' num2str(mean(y.*data_confidences < 0)) '   w: ' num2str(max(w))]);  
        
        classifier.wcs(t, c).dt = dt;
        classifier.wcs(t, c).confidences = confidences;       
             
        
        aveconf(t) = mean(1 ./ (1+exp(-y.*data_confidences)));
        if t>10 && (aveconf(t)-aveconf(t-10) < stopval)
            disp(num2str(aveconf))
            disp(['Stopping after ' num2str(t) ' trees'])            
            break;
        end
        
    end

    finalconf = 1 ./ (1+exp(-y.*data_confidences));
    finalerr = (y.*data_confidences < 0);
    disp(['confidence:: mean: ' num2str(mean(finalconf)) ...
        '  pos: ' num2str(mean(finalconf(y==1))) ...
        '  neg: ' num2str(mean(finalconf(y~=1)))]);
    disp(['training error:: mean: ' num2str(mean(finalerr)) ...
        '  pos: ' num2str(mean(finalerr(y==1))) ...
        '  neg: ' num2str(mean(finalerr(y~=1)))]);    
    all_conf(:, c) = data_confidences+classifier.h0(c);
  
end

% compute and display training error
[tmp, assigned_label] = max(all_conf, [], 2);
conf_matrix = zeros(num_classes, num_classes);
for c = 1:num_classes    
    indices = find(gid==c);
    for c2 = 1:num_classes
        conf_matrix(c, c2) = mean(assigned_label(indices)==c2);
    end
    disp([gn{c} ' error: ' num2str(mean(assigned_label(indices)~=c))]);
end
disp('Confusion Matrix: ');
disp(num2str(conf_matrix));
disp(['total error: ' num2str(mean(assigned_label~=gid))]);


        