function [err, commonSize, numElements] = statsizechk(nparams,varargin)
%STATSIZECHK Check for compatible array sizes.
%   [ERR,COMMONSIZE,NUMELEMENTS] = STATSIZECHK(NPARAMS,A,B,...,M,N,...) or
%   [ERR,COMMONSIZE,NUMELEMENTS] = STATSIZECHK(NPARAMS,A,B,...,[M,N,...])
%   in effect computes size( A + B + ... + zeros(M,N,...) ), and catches
%   any size mismatches.  NPARAMS is the number of array input arguments.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/01/24 09:36:35 $
%
%   Mex file.

% try
%     tmp = 0;
%     for argnum = 1:nparams
%         tmp = tmp + varargin{argnum};
%     end
%     if nargin > nparams+1
%         tmp = tmp + zeros(varargin{nparams+1:end});
%     end
%     err = 0;
%     commonSize = size(tmp);
%     numElements = numel(tmp);
%
% catch
%     err = 1;
%     commonSize = [];
%     numElements = 0;
% end
