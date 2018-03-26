function [doog_filters, texton_data] = APPgetImageFilters;
% Gets filters for creating texture features

% create pair data
angle_step = 15;
filter_size = 5;
sigma = 1.5;
r = 0.25;

% create difference of oriented gaussian filters
num_angles=180/angle_step;
doog_filters = zeros(filter_size, filter_size, num_angles);
for (theta=0:angle_step:180-angle_step)
    doog_filters(:, :, theta/angle_step+1)=createDoogFilter(sigma, r, theta, filter_size); 
end

% load textons
if 0 % DWH no longer uses textons
texton_data = load('unitex_6_1_2_1.4_2_32.mat'); % from Berkeley database
ordered_textons = select_diverse_textons(texton_data.tsim, 12);
new_tims = cell(length(ordered_textons), 1);
for i = 1:length(new_tims)
    new_tims{i} = texton_data.tim{ordered_textons(i)};
end
texton_data.tim = new_tims;
end

texton_data.tim = [];
%num_textons = length(texton_data.tim);