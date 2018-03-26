function [news,errid,errmsg,newrows]=dfgetuserdists(olds,userfun)
%GETUSERDISTS Get user-defined distributions for dfittool
%   [NEWS,ERRID,ERRMSG,NEWROWS]=GETUSERDISTS(OLDS) appends user-defined
%   distribution information to the existing distribution information
%   in the structure OLDS and returns the combined information in the
%   structure NEWS.  Any error id or message is also returned.  NEWROWS
%   is a vector of the indices that are new or have changed.
%
%   [...]=GETUERDISTS(OLDS,USERFUN) uses the function USERFUN in
%   place of the default @dfittooldists.

%   $Revision: 1.1.6.3 $  $Date: 2004/03/02 21:49:25 $
%   Copyright 2003-2004 The MathWorks, Inc.

% If no user-defined function given, use the default if it exists
if nargin<2
   if exist('dfittooldists','file') 
      userfun = @dfittooldists;
   else
      userfun = [];
   end
end

errid = '';
errmsg = '';
news = olds;
newrows = [];

if isempty(userfun)
   return
end

% First try running the user's function
if isa(userfun,'function_handle')
   userfunname = func2str(userfun);
else
   userfunname = userfun;
end
try
   s = feval(userfun);
catch
   errid = 'stats:dfittool:BadUserDistributions';
   errmsg = sprintf(...
       'Error running %s to get user-defined distributions:\n%s',...
                    userfunname, lasterr);
   return
end
if ~isempty(errmsg) || isempty(s)
   return;
end

% Next make sure the result is a structure
if ~isstruct(s)
   errid = 'stats:dfittool:StructureRequired';
   errmsg = sprintf('%s not return a structure',userfunname);;
   return;
end

newfields = fieldnames(s);
numnewdists = length(s);
numolddists = length(olds);
requiredfields = {'name' 'pnames' 'cdffunc' 'pdffunc' 'invfunc'};

% Next make sure the result has all required fields
for j=1:length(requiredfields)
   if isempty(strmatch(requiredfields{j},newfields,'exact'))
      errid = 'stats:dfittool:MissingField';
      errmsg = sprintf('Missing field ''%s'' in %s structure',...
                       requiredfields{j},userfunname);
      return
   end
end

% Make sure the field values are as expected
checked = cell(1,numnewdists);
for j=1:numnewdists
   sj = s(j);
   [errid,errmsg,sj] = checkfields(sj);
   if ~isempty(errid)
      return
   end
   checked{j} = sj;
end

% See if we are overwriting existing fields
if numolddists>0
   oldnames = {olds.name};
   oldcodes = {olds.code};
else
   oldnames = cell(0);
   oldcodes = cell(0);
end

for j=1:numnewdists
   % See if the proposed name or code is in use already
   sj = checked{j};
   newfields = fieldnames(sj);  % may need updating since previous assignment
   name = sj.name;
   code = sj.code;
   oldnrow = strmatch(name,oldnames,'exact');
   oldcrow = strmatch(code,oldcodes,'exact');

   if isempty(oldcrow)
      if ~isempty(oldnrow)
         newrow = oldnrow;          % replace distribution with same name
      else
         newrow = numolddists+1;    % new distribution
      end
   else
      % Trying to re-define an existing distribution
      if ~isempty(oldnrow) && ~isequal(oldcrow,oldnrow)
         errid = 'stats:dfittool:DuplicateName';
         errmsg = sprintf(...
             ['Distribution with code ''%s'' has a name duplicating that ' ...
              'of another distribution.'],code);
         return
      end
      newrow = oldcrow;
   end

   % Update fields in old structure.
   % Can't concatenate with [] if field names differ.
   for fieldnum = 1:length(newfields)
      fieldname = newfields{fieldnum};
      olds(newrow).(fieldname) = sj.(fieldname);
   end

   % Update arrays to guard against duplicates within the new structure
   if newrow>numolddists
      oldnames = [oldnames {name}];
      oldcodes = [oldcodes {code}];
      numolddists = numolddists+1;
   end
   newrows = [newrows newrow];
end

% Return updated structure as new structure
news = olds;
      
% ------------------------------------------
function [errid,errmsg,sj] = checkfields(sj)
%CHECKFIELDS Check that a distribution structure's fields are all valid

% Check required fields
testnames = {'name'};
for j=1:length(testnames)
   field = testnames{j};
   [errid,errmsg,val] = checkstring(sj,field,'',false);
   if ~isempty(errmsg)
      return
   end
   sj.(field) = val;
end

field = 'pnames';
[errid,errmsg,val] = checktext(sj,field,'',false,[]);
if ~isempty(errmsg)
   return
end
sj.(field) = val;
nparams = length(val);

testnames = {'pdffunc' 'cdffunc' 'invfunc'};
for j=1:length(testnames)
   field = testnames{j};
   [errid,errmsg,val] = checkfunc(sj,field,false);
   if ~isempty(errmsg)
      return
   end
   sj.(field) = val;
end

% Check optional fields and fill in defaults
testnames = {'code'};
defaults  = {lower(sj.name)};
for j=1:length(testnames)
   field = testnames{j};
   [errid,errmsg,val] = checkstring(sj,field,defaults{j},true);
   if ~isempty(errmsg)
      return
   end
   sj.(field) = val;
end

testnames = {'hasconfbounds' 'iscontinuous' 'islocscale' 'uselogpp' ...
             'censoring'     'paramvec'};
defaults  = {false           true           false        false      ...
             false           true};
for j=1:length(testnames)
   field = testnames{j};
   default = defaults{j};
   [errid,errmsg,val] = checklogical(sj,field,default);
   if ~isempty(errmsg)
      return
   end
   sj.(field) = val;
end

testnames = {'likefunc' 'logcdffunc' 'loginvfunc'};
for j=1:length(testnames)
   field = testnames{j};
   [errid,errmsg,val] = checkfunc(sj,field,true);
   if ~isempty(errmsg)
      return
   end
   sj.(field) = val;
end

field = 'prequired';
[errid,errmsg,val] = checklogical(sj,field,false(1,nparams),nparams);
if ~isempty(errmsg)
   return
end
sj.(field) = val;

field = 'pdescription';
[errid,errmsg,val] = checktext(sj,field,{},true,nparams);
if ~isempty(errmsg)
   return
end
sj.(field) = val;

field = 'closedbound';
[errid,errmsg,val] = checklogical(sj,field,[false false],2);
if ~isempty(errmsg)
   return
end
sj.(field) = val;

field = 'support';
[errid,errmsg,val] = checksupport(sj);
if ~isempty(errmsg)
   return
end
sj.(field) = val;


% ------------------------------------------
function [errid,errmsg,val] = checklogical(s,field,default,nvals)
%CHECKLOGICAL Check that a field has a valid logical value

if nargin<4
   nvals = 1;
end

val = [];
errid = '';
errmsg = '';
if ~isfield(s,field) || isempty(s.(field))
   val = default;
else
   val = s.(field);
end

if (numel(val) ~= nvals)
   errid = 'stats:dfittool:WrongSize';
   errmsg = sprintf(...
       'The ''%s'' field must contain %d element(s).',field,nvals);
elseif ~(islogical(val) || isnumeric(val))
   errid = 'stats:dfittool:NotLogical';
   errmsg = sprintf(...
       'The ''%s'' field must be true or false.',field);
else
   val = (val ~= 0);
end


% ------------------------------------------
function [errid,errmsg,val]=checktext(s,field,default,optional,nvals)
%CHECKTEXT Check that a field has a value that is an array of strings

if nargin<5
   nvals = 1;
end

errid = '';
errmsg = '';
val = '';
if ~isfield(s,field) || isempty(s.(field))
   if optional
      val = default;
      return
   else
      errid = 'stats:dfittool:EmptyNotAllowed';
      errmsg = sprintf('The ''%s'' field must not be empty.',field);
      return
   end
end
val = s.(field);
if iscellstr(val)
   if ~isempty(nvals) && numel(val)~=nvals
      errid = 'stats:dfittool:BadSize';
      errmsg = sprintf(...
          'The ''%s'' field must contain %d string(s).',field,nvals);
   end
elseif ischar(val)
   if ~isempty(nvals) && size(val,1)~=nvals
      errid = 'stats:dfittool:BadSize';
      errmsg = sprintf(...
          'The ''%s'' field must contain %d string(s).',field,nvals);
   else
      val = cellstr(val);
   end
else
   errid = 'stats:dfittool:NotCharacter';
   errmsg = sprintf('Value in ''%s'' field must be a character array or cell array of strings.',field);
end

% ------------------------------------------
function [errid,errmsg,val]=checkstring(s,field,default,optional)
%CHECKSTRING Check that a field has a valid string

errid = '';
errmsg = '';
val = '';
if ~isfield(s,field) || isempty(s.(field))
   if optional
      val = default;
      return
   else
      errid = 'stats:dfittool:EmptyNotAllowed';
      errmsg = sprintf('The ''%s'' field must not be empty.',field);
      return
   end
end
val = s.(field);
if ~ischar(val) || (~isequal(size(val), [1,length(val)]))
   errid = 'stats:dfittool:NotCharacter';
   errmsg = sprintf('Value in ''%s'' field must be a character string.',field);
end


% ------------------------------------------
function [errid,errmsg,val] = checkfunc(s,field,optional)
%CHECKFUNC Check that a field has a valid function value


val = '';
errid = '';
errmsg = '';
if ~isfield(s,field) || isempty(s.(field))
   if ~optional
      errid = 'stats:dfittool:EmptyNotAllowed';
      errmsg = sprintf('The ''%s'' field must not be empty.',field);
   end
   return
end

val = s.(field);
if ~isa(val,'function_handle') || ~isscalar(val)
   errid = 'stats:dfittool:NotFunctionHandle';
   errmsg = sprintf(...
          'The ''%s'' field must contain a single function handle.',field);
end

% ------------------------------------------
function [errid,errmsg,val] = checksupport(s,default)
%CHECKLOGICAL Check that a field has a valid logical value

val = [];
errid = '';
errmsg = '';
default = [-Inf Inf];
field = 'support';
if ~isfield(s,field) || isempty(s.(field))
   val = default;
   return;
end
val = s.(field);
if ~isnumeric(val) || numel(val)~=2
   errid = 'stats:dfittool:BadSupport';
   errmsg = 'The ''support'' value must contain a two-element vector.';
elseif (val(1)>=val(2)) || any(isnan(val))
   errid = 'stats:dfittool:BadSupport';
   errmsg = 'The ''support'' values must be increasing and non-NaN.';
end
