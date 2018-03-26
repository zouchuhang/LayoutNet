function feature = HOGMe(image, patchSize, gridSpacing)

% if no input paramaters, return the dimension of the feature
if nargin == 0
  feature = 31;
  return
end

if patchSize == gridSpacing
    % just standard Pedro's HOG
    feature = features_pedro(image,patchSize);
else
    % space and grid are different
end