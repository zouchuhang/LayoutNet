function varargout = plot3D(varargin)
% downloaded from http://www.mathworks.com/matlabcentral/fileexchange/21044-3d-voxelizer
% PLOT3D Draws a collection of voxel points in 3 -Dimensional space.
%    PLOT3D(A) Draws the voxels contained in A. A is a 3-dimensional matrix
%       of zeros and ones, where 1-valued elements of the matrix will be
%       drawn and 0-elements will be left blank. The (I,J,K) element of the
%       matrix will be drawn at the coordinates (I,J,K). Alternatively, A
%       can be a Mx3 vector of vertex points [ x y z ].
%    PLOT3D(A,SIZE) as above. SIZE defines the drawing voxels edge length.
%    PLOT3D(A,SIZE,C) as above. C defines the drawing colour of the voxels.
%    PLOT3D(A,SIZE,C,STYLE) as above. STYLE defines the style of marker to be
%       used.
% 
%       Possible marker styles are:
%          .   point
%          o   circle
%          x   x-mark
%          *   star
%          s   square
%          d   diamond
%          v   triangle (down)
%          ^   triangle (up)
%          <   triangle (left)
%          >   triangle (right)
%          p   pentagram
%          h   hexagram
%          vox voxel
% 
%    PLOT3D(A,SIZE,C,STYLE,ALPHA) as above. ALPHA defines the transparenciy
%       of the voxels.
%    PLOT3D(...,MODE) MODE defines the interactive drawing mode to be used.
%       This is intended to aid visualization. Available modes are: 
% 
%          pasive    Draws all in one go
%          timed     Draws next interval every X seconds. Time interval X
%                    is set is passed as an additional parameter.
%          keyboard  Draws next layer when keyboard is hit
%
%    h = PLOT3D(...) Returns a vector containing the handles of all drawn 
%       voxels.
% 
%    Example 1
%    ---------
%    Creates the occupancy matrix of a discretized sphere with radius 10.
% 
%       r = 10;
%       side = 30;
%       sphere = zeros(side,side,side);
%       for i=1:side; x = i-side/2;
%         for j=1:side; y = j-side/2;
%           for k=1:side; z = k-side/2;
%             if sum([ x y z ].^2)<=r^2; sphere(i,j,k) = 1; end;
%           end
%         end
%       end
%       plot3D(sphere,'pasive');
%       plot3D(Sphere_1,'timed', 0.1);
%       plot3D(Sphere_1,'keyboard');
%
% See also VOXEL

  % Constants
  MODE_1 = 'pasive';
  MODE_2 = 'timed';
  MODE_3 = 'keyboard';

  c2 = 'c';
  timed_t = 3;
  VOXEL = 'vox';
  axis vis3d;
  view([30 30]);

  % Default values
  vox_size = 1;
  c = 'm';
  alpha = 1;
  style = VOXEL;
  mode = MODE_2;

  if nargin>=1
    tmp_mode = varargin{length(varargin)};
    tmp_t = NaN;
    if nargin>=2 && isnumeric(tmp_mode)
      tmp_t = tmp_mode;
      tmp_mode = varargin{length(varargin)-1};
    end

    if strcmpi(tmp_mode,MODE_1) || strcmpi(tmp_mode,MODE_2) || ...
            strcmpi(tmp_mode,MODE_3)
      varargin(length(varargin)) = [];
      mode = tmp_mode;

      if ~isnan(tmp_t)
        varargin(length(varargin)) = [];
        timed_t = tmp_t;
      end
    end
  end

  if isempty(varargin)
    disp('Self test!!');
    a = zeros(10,10,10);
    for i=1:11; x = i-6;
      for j=1:11; y = j-6;
        for k=1:11; z = k-6;
          if sum([ x y z ].^2)<=25; a(i,j,k) = 1; end;
        end
      end
    end
  end

  if length(varargin)>=1;
    a = varargin{1};
  end;
  if length(varargin)>=2;
    vox_size = varargin{2};
  end;
  if length(varargin)>=3;
    c = varargin{3};
  end;
  if length(varargin)>=4;
    style = varargin{4};
  end;
  if length(varargin)>=5;
    alpha = varargin{5};
  end;


  if ndims(a)==2
    % Already in coordinate format
    if strcmp(mode,MODE_2) || strcmp(mode,MODE_3)
      S = a;
    else
      S = [];
    end
    S2 = S;  % optimization required for coord format case
  elseif ndims(a)==3
    % 3D occupancy matrix --> convert to coordinate format
    if strcmp(mode,MODE_2) || strcmp(mode,MODE_3)
      S = sparse3D(a);
    else
      S = [];
    end

    if license('test','image_toolbox')
      b = bwperim(a);
      S2 = sparse3D(b);
    else
      S2 = sparse3D(a);
    end
  else
    error('Missmatching numer of dimensions: must be 2 or 3');
  end

  % Start drawing
  box = [ 1 1 1 ] * vox_size;

  hold on;
  h = [];
  zs = unique(S2(:,3));
  for i=1:length(zs)

    % Pre-draw
    if strcmp(mode,MODE_2) || strcmp(mode,MODE_3)
      idx  = find(S(:,3)==zs(i));
      hd1 = [];
      if strcmpi(style, VOXEL)
        for j=1:size(idx,1)
          hd2 = voxel(S(idx(j),:),box,c2,alpha);
          hd1 = vertcat(hd1,hd2);
        end
      else
        hd1 = plot3(S(idx,1), S(idx,2), S(idx,3), [c2 style]);
      end

      if strcmp(mode,MODE_2)
        pause(timed_t);
      else
        pause;
      end
      delete(hd1);
    end

    idx2 = find(S2(:,3)==zs(i));
    h1 = [];
    if strcmpi(style, VOXEL)
      for j=1:size(idx2,1)
        h2 = voxel(S2(idx2(j),:),box,c,alpha);
        h1 = vertcat(h1,h2);
      end
    else
      h1 = plot3(S2(idx,1), S2(idx,2), S2(idx,3), [c style]);
    end
    h = vertcat(h1,h);
  end

  if nargout>0
    varargout{1} = h;
  end