function name = dfgetfitname()

% $Revision: 1.1.6.2 $  $Date: 2004/01/24 09:35:38 $
% Copyright 2003-2004 The MathWorks, Inc.

count=dfgetset('fitcount');
if isempty(count)
    count = 1;
end
taken = 1;
while taken
    name=sprintf('fit %i', count);
    if isempty(find(getfitdb,'name',name))
        taken = 0;
    else
        count=count+1;
    end
end
dfgetset('fitcount',count+1);

