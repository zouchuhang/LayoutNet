function a = on2off(b)
%ON2OFF Simple helper returns 'on' given 'off' and vice versa

%   $Revision: 1.1.6.2 $  $Date: 2004/01/24 09:36:22 $
%   Copyright 2003-2004 The MathWorks, Inc.

if isequal(b,'on')
   a = 'off';
else
   a = 'on';
end
