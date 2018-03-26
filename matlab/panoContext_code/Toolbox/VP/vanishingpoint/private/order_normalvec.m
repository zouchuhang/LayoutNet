function normal_vec = order_normalvec(normal_vec)
% Roughly orders normal_vec.
% This is just approximate as normal_vec's may be aligned with the camera orientation.
% 1st normal_vec: the one orthogonal to the ground plane (up-down)
% 2nd normal_vec: the one going left-right
% 3rd normal_vec: the one parallel to the direction that the camera is looking at (Z direction)
% In other words, 1st: Y, 2nd: X, 3rd: Z, where XYZ is
%                  > Z (camera direction)
%                 /
%                /
%               --------> X
%               |
%               |
%               |
%               v Y
% 

n1 = normal_vec{1}; % assume vanishline.m is correct.
n2 = normal_vec{2};
n3 = normal_vec{3};

if abs([1 0 0] * n2) > abs([1 0 0] * n3)
    normal_vec{2} = n2;
    normal_vec{3} = n3;
else
	% swap
    normal_vec{2} = n3;
    normal_vec{3} = n2;
end

