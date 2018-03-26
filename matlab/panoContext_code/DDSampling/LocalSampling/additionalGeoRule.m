function [ score ] = additionalGeoRule( obj_xyz, score, val )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

floatValid = obj_xyz(:,3)>-150;
underValid = all(obj_xyz(:,1:2)<0,2) & all(obj_xyz(:,4:5)>0,2);

%% bed(9), chair(10), sofa(12) cannot float
score(floatValid,[9 10 12]) = val;

%% bedside(2), desk(7), wardrobe(11), cabinet(5)
score(floatValid,[5]) = val;

%% nothing under foot
score(underValid,:) = val;

end

