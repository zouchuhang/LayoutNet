function features = APPgetVpFeatures(spinfo, lines, rCenter, imsize)
% Get the vanishing point features for superpixels defined by spinfo
% spinfo(nSp).{lines(nLinesInSp)}
% rCenter([y x]): region center wrt image center (-0.5 to 0.5)
%
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  09/30/2005

features = zeros(1, 19); 

num_blocks = length(spinfo);
used_lines = zeros(size(lines, 1), 1);
max_angle = pi/8;

for i = 1:num_blocks
    used_lines(spinfo(i).lines) = 1;
end
lines = lines(find(used_lines>0), :);
num_lines = size(lines, 1);

features(:) = APPlines2vpFeatures(lines, max_angle);
