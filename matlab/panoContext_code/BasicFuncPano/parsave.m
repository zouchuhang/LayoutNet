function parsave( filename, varargin  )
%PARSAVE Summary of this function goes here
%   Detailed explanation goes here
savestring = [];
for i = 1:2:length(varargin)
    eval(sprintf('%s=varargin{%d};', varargin{i}, i+1));
    savestring = [savestring ',''' varargin{i} ''''];
end
eval(['save( ''' filename ''' ' savestring ');']);
end

