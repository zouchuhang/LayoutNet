function [new, fittype]=dfcreatecopy(original);

%   $Revision: 1.1.6.3 $  $Date: 2004/01/24 09:35:23 $
%   Copyright 2003-2004 The MathWorks, Inc.

fittype = original.fittype;

new = copyfit(original);
new = java(new);

