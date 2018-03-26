%% pattern
num_ap = 25;

ap = exp(linspace(log(0.2), log(5), num_ap)); %height/width
width = zeros(num_ap*2, 1);
height = zeros(num_ap*2, 1);
for i = 0:num_ap-1
    if ap(i+1)>1
        width(i*2+1) = 10;
        width(i*2+2) = 12;
        height(i*2+1) = 10*ap(i+1);
        height(i*2+2) = 12*ap(i+1);
    else
        height(i*2+1) = 10;
        height(i*2+2) = 12;
        width(i*2+1) = 10/ap(i+1);
        width(i*2+2) = 12/ap(i+1);
    end
end
pattern = zeros(num_ap*2, 18);
patternSz = zeros(num_ap*2, 2);
for i = 1:num_ap*2
    w = round(width(i));   h = round(height(i));
    ws = [1 round(w/2)+1 w+1];
    hs = [1 round(h/2)+1 h+1];
    pattern(i,1:2:end) = ws([1 2 3 1 2 3 1 2 3]);
    pattern(i,2:2:end) = hs([1 1 1 2 2 2 3 3 3]);
    patternSz(i,:) = [max(pattern(i,1:2:end))+4 max(pattern(i,2:2:end))+4];
end

%%
config.pattern = pattern;
config.patternSz = patternSz;
config.nmsThresh = 0.75;
config.repThresh = -0.25;
config.max_scale = 1;
config.min_scale = 0.02;
config.newRepThresh = -0.25;

%% detection file
config.modelfile = './ObjectHypothesisGeneration/RectangleDetector/finalModel_4.mat';