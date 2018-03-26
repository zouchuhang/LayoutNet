function feature = feature_pyramid(image, detect_max_scale, detect_min_scale)

%Make sure image is in double format
image = double(image);
if nargin==1
    detect_max_scale = 2.0;
    detect_min_scale = .25;
% detect_max_scale = 2.0;
% detect_min_scale = 0.1;
end

%Hardcoded maximum number of levels in the pyramid
MAXLEVELS = 200;

%Hardcoded minimum dimension of smallest (coarsest) pyramid level
MINDIMENSION = 5;

%Get the levels per octave from the parameters
interval = 5;

sc = 2 ^(1/interval);

% Start at detect_max_scale, and keep going down by the increment sc, until
% we reach MAXLEVELS or detect_min_scale
feature.scale = zeros(1,MAXLEVELS);
feature.feat = {};
for i = 1:MAXLEVELS
    scaler = detect_max_scale / sc^(i-1);
    
    if scaler < detect_min_scale
        return
    end
    
    feature.scale(i) = scaler;
    if feature.scale(i)>1
        scaled = imresize(image, feature.scale(i));
    else
        scaled = resize(image,feature.scale(i));
    end
    %if minimum dimensions is less than or equal to 5, exit
    if min([size(scaled,1) size(scaled,2)])<=MINDIMENSION
        feature.scale = feature.scale(feature.scale>0);
        return;
    end
    
    feature.feat{i} = HOGMe(scaled,8,8);
    
    %if we get zero size feature, backtrack one, and dont produce any
    %more levels
    if (size(feature.feat{i},1)*size(feature.feat{i},2)) == 0
        feature.feat = feature.feat(1:end-1);
        feature.scale = feature.scale(1:end-1);
        return;
    end
    
    %recover lost bin!!!
    feature.feat{i} = padarray(feature.feat{i}, [1 1 0], 0);
    
    %if the max dimensions is less than or equal to 5, dont produce
    %any more levels
    if max([size(feature.feat{i},1) size(feature.feat{i},2)])<=MINDIMENSION
        feature.scale = feature.scale(feature.scale>0);
        return;
    end
end