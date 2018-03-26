function [emsg,link,dlink,ilink,exponent]=stattestlink(link)
%STATTESTLINK  Test link function for GLMFIT and GLMVAL

%   Author:  Tom Lane, 3-7-2000
%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 1.3 $  $Date: 2002/02/04 19:25:49 $

dlink = '';
ilink = '';
exponent = {};
emsg = '';

if (iscell(link))
   % A cell array of three functions is okay
   if (length(link)~=3)
      emsg = 'LINK cell array must have three components';
      return
   end
   dlink = link{2};
   ilink = link{3};
   link = link{1};
   if (~statglmeval('testlink',link))
      emsg = 'LINK function is not valid';
      return
   end
   if (~statglmeval('testlink',dlink))
      emsg = 'LINK function derivative is not valid';
      return
   end
   if (~statglmeval('testlink',ilink))
      emsg = 'LINK function inverse is not valid';
      return
   end
elseif (ischar(link) & size(link,1)==1)
   % A function name is okay, but three functions must exist
   dlink = ['d_' link];
   ilink = ['i_' link];
   if (~statglmeval('testlink',link))
      emsg = sprintf('Cannot find LINK function %s.',link);
      return
   end
   if (~statglmeval('testlink',dlink))
      emsg = sprintf('Cannot find LINK function derivative %s.',dlink);
      return
   end
   if (~statglmeval('testlink',ilink))
      emsg = sprintf('Cannot find LINK function inverse %s.',ilink);
      return
   end
elseif (isnumeric(link) & length(link)==1)
   exponent = {link};
   link = 'power';
   dlink = ['d_' link];
   ilink = ['i_' link];
else
   emsg = 'LINK function is not valid.';
   return
end
