function [emsg,varargout]=statgetargsuser(pnames,dflts,delim,varargin)
%STATGETARGS Process parameter name/value pairs for statistics functions
%   [EMSG,A,B,... USER]=
%       STATGETARGS(PNAMES,DFLTS,DELIM,'NAME1',VAL1,'NAME2',VAL2,...)
%   accepts a cell array PNAMES of valid parameter names, a cell array
%   DFLTS of default values for the parameters named in PNAMES, the
%   name of a "user parameter delimiter", and additional parameter name/value
%   pairs.  Returns parameter values A,B,... in the same order as the names
%   in PNAMES, plus a single cell array USER of all input args that follow
%   the name specified by DELIM.  Outputs corresponding to entries in PNAMES
%   that are not specified in the name/value pairs are set to the
%   corresponding value from DFLTS, and USER if set to [] if the name in
%   DELIM does not appear in the name/value pairs.  If nargout is equal to
%   length(PNAMES)+2, then unrecognized name/value pairs are an error.  If
%   nargout is equal to length(PNAMES)+3, then all unrecognized name/value
%   pairs are returned in a single cell array following any other outputs.
%
%   EMSG is empty if the arguments are valid, or the text of an error message
%   if an error occurs.  STATGETARGS does not actually throw any errors, but
%   rather returns an error message so that the caller may throw the error.
%   Outputs will be partially processed after an error occurs.
%
%   This utility is used by some Statistics Toolbox functions to process
%   name/value pair arguments.
%
%   Example:
%       pnames = {'color' 'linestyle', 'linewidth'}
%       dflts  = {    'r'         '_'          '1'}
%       delim = 'userargs'
%       varargin = {'linew' 2 'nonesuch' [1 2 3] 'linestyle' ':' ...
%          'userargs' 'pretty much' 'anything here' [1 2 3] {1 2 3}}
%       [emsg,c,ls,lw,ua] = statgetargs(pnames,dflts,varargin{:})    % error
%       [emsg,c,ls,lw,ua,ur] = statgetargs(pnames,dflts,varargin{:}) % ok
%       % c is 'r', ls is ':', and lw is 2
%       % ua is {'pretty much' 'anything here' [1 2 3] {1 2 3}}
%       % ur is {'nonesuch' [1 2 3]}

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 1.2 $  $Date: 2002/02/04 19:25:45 $

% We always create (nparams+2) outputs:
%    one for emsg
%    nparams varargs for values corresponding to names in pnames
%    one vararg for all user values
% If they ask for one more (nargout == nparams+3), it's for unrecognized
% names/values

% Initialize some variables
emsg = '';
nparams = length(pnames);
varargout = {dflts{:} {}}; % user args default to empty
unrecog = {};
nargs = length(varargin);

pnames = {pnames{:} delim}; % (nparams+1)th element is the user args delimiter

% Process name/value pairs
for j=1:2:nargs
    pname = varargin{j};
    if ~ischar(pname)
        emsg = sprintf('Parameter name must be text.');
        break;
    % Must have one or more args after a name
    elseif j == nargs
        emsg = sprintf('Wrong number of arguments.');
        break;
    end
    i = strmatch(lower(pname),pnames);
    if isempty(i)
        % if they've asked to get back unrecognized names/values, add this
        % one to the list
        if nargout > nparams+2
            unrecog((end+1):(end+2)) = {varargin{j} varargin{j+1}};
            
        % otherwise, it's an error
        else
            emsg = sprintf('Invalid parameter name:  %s.',pname);
            break;
        end
    elseif length(i)>1
        emsg = sprintf('Ambiguous parameter name:  %s.',pname);
        break;
        
    % matched delimiter, everything remaining is user args, return them
    % all untouched as a single cell array
    elseif i == nparams+1
        varargout{i} = varargin((j+1):end);
        break;
    
    % a regular old name/value pair
    else
        varargout{i} = varargin{j+1};
    end
end

varargout{nparams+2} = unrecog;
