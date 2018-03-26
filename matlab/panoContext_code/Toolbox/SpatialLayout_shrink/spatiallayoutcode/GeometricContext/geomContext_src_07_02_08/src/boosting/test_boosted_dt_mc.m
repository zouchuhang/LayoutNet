function confidences = test_boosted_dt_mc(classifier, features)
% confidences = test_boosted_dt_mc(classifier, features)
%
% Returns a log likelihod ratio for each class in the classifier    
% 
% Input:
%  classifier: boosted decision tree classifier
%  features:   classifier features (ndata, nvariables)
% Output:
%   confidences(ndata, nclasses): 
%      P(class=k|features) \propto 1./(1+exp(-confidences(k)))

npred = classifier.wcs(1).dt.npred;
if size(features, 2)~=npred
    error('Incorrect number of attributes')
end

wcs = classifier.wcs;  
nclasses = size(wcs, 2);

ntrees = size(wcs, 1);

confidences = zeros(size(features, 1), nclasses);
for c = 1:nclasses    
    for t = 1:ntrees        
        if ~isempty(wcs(t,c).dt)                                              
            if 1
                dt = wcs(t,c).dt; 
%                 dt = tree_getNewVersion(dt);
                [var, cut, children, catsplit] = tree_getParameters(dt);
                nodes = treevalc(int32(var), cut, int32(children(:, 1)), ...
                        int32(children(:, 2)), catsplit(:, 1), features');  
                %disp(num2str(nodes));      
            else
                [class_indices, nodes, classes] = treeval(wcs(t, c).dt, features);             
            end
            confidences(:, c) = confidences(:, c) + wcs(t, c).confidences(nodes);
        end        
    end
    confidences(:, c) = confidences(:, c) + classifier.h0(c);
end

   
