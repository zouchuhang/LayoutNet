%            3  
%1'----------------------2'
%  |                    |
%  |                    |
% 1|                    |2
%  |                    |
%  |                    |
%4'----------------------3'
%            4

%% it is better to define by mapping cuboid to rectangle
GeoRule(1).viewID = 1;
GeoRule(1).surfaceID = 2;
GeoRule(1).pointID = [6 8 4 2];
GeoRule(2).viewID = 2;
GeoRule(2).surfaceID = 1;
GeoRule(2).pointID = [7 5 1 3];
GeoRule(3).viewID = 3;
GeoRule(3).surfaceID = 4;
GeoRule(3).pointID = [8 7 3 4];
GeoRule(4).viewID = 4;
GeoRule(4).surfaceID = 3;
GeoRule(4).pointID = [5 6 2 1];
GeoRule(5).viewID = 5;
GeoRule(5).surfaceID = 6;
GeoRule(5).pointID = [7 8 6 5];
GeoRule(6).viewID = 6;
GeoRule(6).surfaceID = 5;
GeoRule(6).pointID = [1 2 4 3];

CubRule(1).suf = [1 3];
CubRule(1).cont = [4 8];
CubRule(1).check = 5;

CubRule(2).suf = [4 1];
CubRule(2).cont = [2 6];
CubRule(2).check = 5;

CubRule(3).suf = [2 3];
CubRule(3).cont = [3 7];
CubRule(3).check = 5;

CubRule(4).suf = [2 4];
CubRule(4).cont = [1 5];
CubRule(4).check = 5;

CubRule(5).suf = [5 1];
CubRule(5).cont = [6 8];
CubRule(5).check = [3 4];

CubRule(6).suf = [5 2];
CubRule(6).cont = [5 7];
CubRule(6).check = [3 4];

CubRule(7).suf = [3 5];
CubRule(7).cont = [7 8];
CubRule(7).check = [1 2];

CubRule(8).suf = [4 5];
CubRule(8).cont = [5 6];
CubRule(8).check = [1 2];

%% cuboid rule version 1
CuboidRule(1).lhs.vid = 1;
CuboidRule(1).lhs.lid = 2;
CuboidRule(1).lhs.pid = [3 2];
CuboidRule(1).rhs.vid = 3;
CuboidRule(1).rhs.lid = 1;
CuboidRule(1).lhs.pid = [4 1];
CuboidRule(1).chk.condition = [2 3 0];
CuboidRule(1).chk.vid = 5;
CuboidRule(1).chk.lhslid = [3 2];
CuboidRule(1).chk.rhslid = [3 3];

CuboidRule(2).lhs.vid = 3;
CuboidRule(2).lhs.lid = 2;
CuboidRule(2).rhs.vid = 2;
CuboidRule(2).rhs.lid = 1;
CuboidRule(2).chk.condition = [3 0];
CuboidRule(2).chk.vid = 5;
CuboidRule(2).chk.lhslid = [3 3];
CuboidRule(2).chk.rhslid = [3 1];

CuboidRule(3).lhs.vid = 4;
CuboidRule(3).lhs.lid = 2;
CuboidRule(3).rhs.vid = 1;
CuboidRule(3).rhs.lid = 1;
CuboidRule(3).chk.condition = [3 0];
CuboidRule(3).chk.vid = 5;
CuboidRule(3).chk.lhslid = [3 4];
CuboidRule(3).chk.rhslid = [3 2];

CuboidRule(4).lhs.vid = 2;
CuboidRule(4).lhs.lid = 2;
CuboidRule(4).rhs.vid = 4;
CuboidRule(4).rhs.lid = 1;
CuboidRule(4).chk.condition = [3 0];
CuboidRule(4).chk.vid = 5;
CuboidRule(4).chk.lhslid = [3 1];
CuboidRule(4).chk.rhslid = [3 4];

%-------------------------------------
CuboidRule(5).lhs.vid = 1;
CuboidRule(5).lhs.lid = 3;
CuboidRule(5).rhs.vid = 5;
CuboidRule(5).rhs.lid = 2;
CuboidRule(5).chk(1).condition = [2 1];
CuboidRule(5).chk(1).vid = 4;
CuboidRule(5).chk(1).lhslid = [1 2];
CuboidRule(5).chk(1).rhslid = [4 3];
CuboidRule(5).chk(2).condition = [2 0];
CuboidRule(5).chk(2).vid = 3;
CuboidRule(5).chk(2).lhslid = [2 1];
CuboidRule(5).chk(2).rhslid = [3 3];

CuboidRule(6).lhs.vid = 2;
CuboidRule(6).lhs.lid = 3;
CuboidRule(6).rhs.vid = 5;
CuboidRule(6).rhs.lid = 1;
CuboidRule(6).chk(1).condition = [2 1];
CuboidRule(6).chk(1).vid = 4;
CuboidRule(6).chk(1).lhslid = [2 1];
CuboidRule(6).chk(1).rhslid = [3 1];
CuboidRule(6).chk(2).condition = [2 0];
CuboidRule(6).chk(2).vid = 3;
CuboidRule(6).chk(2).lhslid = [1 2];
CuboidRule(6).chk(2).rhslid = [3 3];

CuboidRule(7).lhs.vid = 3;
CuboidRule(7).lhs.lid = 3;
CuboidRule(7).rhs.vid = 5;
CuboidRule(7).rhs.lid = 3;
CuboidRule(7).chk(1).condition = [1 0];
CuboidRule(7).chk(1).vid = 1;
CuboidRule(7).chk(1).lhslid = [1 2];
CuboidRule(7).chk(1).rhslid = [2 3];
CuboidRule(7).chk(2).condition = [1 1];
CuboidRule(7).chk(2).vid = 2;
CuboidRule(7).chk(2).lhslid = [2 1];
CuboidRule(7).chk(2).rhslid = [1 3];

CuboidRule(8).lhs.vid = 4;
CuboidRule(8).lhs.lid = 3;
CuboidRule(8).rhs.vid = 5;
CuboidRule(8).rhs.lid = 4;
CuboidRule(8).chk(1).condition = [1 0];
CuboidRule(8).chk(1).vid = 1;
CuboidRule(8).chk(1).lhslid = [2 1];
CuboidRule(8).chk(1).rhslid = [2 3];
CuboidRule(8).chk(2).condition = [1 1];
CuboidRule(8).chk(2).vid = 2;
CuboidRule(8).chk(2).lhslid = [1 2];
CuboidRule(8).chk(2).rhslid = [1 3];
