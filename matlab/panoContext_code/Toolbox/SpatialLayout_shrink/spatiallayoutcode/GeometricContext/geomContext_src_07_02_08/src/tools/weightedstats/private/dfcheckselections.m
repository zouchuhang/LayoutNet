function [err,d,c,f]=dfcheckselections(data,censoring,frequency,dval,cval,fval)

% For use by DFITTOOL

%   $Revision: 1.1.6.6 $  $Date: 2004/01/24 09:35:21 $
%   Copyright 2003-2004 The MathWorks, Inc.

err = '';
d_l = 0;
c_l = 0;
f_l = 0;
NONE = '(none)';
d = [];
c = [];
f = [];

if isequal(data, NONE)
    err = sprintf('Invalid Data Choice: %s\n', NONE);
else
    if isempty(data)
        dataname = 'data variable';
    else
        dataname = sprintf('"%s"',data);
    end
    try
        if nargin<4
            d=evalin('base',data);
        else
            d = dval;
        end
        if isvector(d) && (length(d) > 1)
            if any(isinf(d))
               err = sprintf('%s cannot contain Inf or -Inf\n', dataname);
            elseif ~isreal(d)
               err = sprintf('%s cannot be complex\n', dataname);
            else
               d_l = length(d);
            end
        else
            err = sprintf('%s is not a vector\n', dataname);
        end
    catch
        err = [err sprintf('Invalid expression: %s\n %s\n', dataname, lasterr)];
    end
end

if (nargin<5) && (isempty(censoring) || isequal(censoring, NONE))
    c_l = -1;
else
    if isempty(censoring)
        censname = 'censoring variable';
    else
        censname = sprintf('"%s"',censoring);
    end
    try
        if nargin<5
            c=evalin('base',censoring);
        else
            c = cval;
        end
        if isempty(c)
           c_l = -1;
        elseif isvector(c) && (length(c) > 1)
            if ~all(ismember(c, 0:1))
                err = [err sprintf('%s must be a logical vector.\n',censname)];
            elseif any(isinf(c))
                err = [err sprintf('%s cannot contain Inf or -Inf\n', censname)];
            elseif ~isreal(c)
                err = [err sprintf('%s cannot be complex\n', censname)];
            else
                c_l = length(c);
            end
        else
            err = [err sprintf('%s is not a vector\n', censname)];
        end
    catch
        err = [err sprintf('Invalid expression: %s\n %s\n', censname, lasterr)];
    end
end

if (nargin<6) && (isempty(frequency) || isequal(frequency, NONE))
    f_l = -1;
else
    if isempty(frequency)
        freqname = 'frequency variable';
    else
        freqname = sprintf('"%s"',frequency);
    end
    try
        if nargin<6
            f=evalin('base',frequency);
        else
            f = fval;
        end
        if isempty(f)
           f_l = -1;
        elseif isvector(f) && (length(f) > 1)
            if any(f<0) || any(f~=round(f) & ~isnan(f))
               err = [err sprintf('%s values must be non-negative integers.\n',freqname)];
            elseif any(isinf(f))
               err = [err sprintf('%s cannot contain Inf or -Inf\n', freqname)];
            elseif ~isreal(f)
               err = [err sprintf('%s cannot be complex\n', freqname)];
            else
               f_l = length(f);
            end
        else
            err = [err sprintf('%s is not a vector\n', freqname)];
        end
    catch
        err = [err sprintf('Invalid expression: %s\n %s\n', freqname, lasterr)];
    end
end

% Check lengths if no other errors
if isequal(err, '')
    if ((c_l ~= -1) && (c_l ~= d_l)) || ((f_l ~= -1) && (f_l ~= d_l))
        err = sprintf('Vector lengths must be equal\n');
        err = [err sprintf('    Data length: %d\n', d_l)];
        if (c_l ~= -1) && (c_l ~= d_l)
            err = [err sprintf('    Censoring length: %d\n', c_l)];
        end
        if (f_l ~= -1) && (f_l ~= d_l)
            err = [err sprintf('    Frequency length: %d\n', f_l)];
        end
    end
end

% Must have some non-censored data
if isempty(err) && c_l~=-1
   if (f_l==-1 && all(c==1)) || (f_l~=-1 && all(c(f>0)==1))
      err = 'Cannot have all observations censored';
   end
end

      