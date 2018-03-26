function P = slparseprops(P0, varargin)
%SLPARSEPROPS Parses input parameters
%
% $ Syntax $
%   - P = slparseprops(P0, property_struct)
%   - P = slparseprops(P0, property_name1, property_value1, ...)
%
% $ Syntax $
%   - P = slparseprops(P0, property_struct) parses the properties from a 
%     property structure. The structure has multiple entries, with the 
%     property names as field names, and the property values are corresponding
%     values. 
%
%   - P = slparseprops(P0, property_name1, property_value1, ...) parses
%     the properties from property list. The list are a argument list with
%     the form of multiple pairs of property names and corresponding values.
%
% $ Remarks $
%   - The P0 specifies the default values of properties. All properties 
%     that can be set should appear in P0. If some properties specified in
%     following arguments are not included in P0, an error will be raised.
%
%   - Passing arguments in the way of property specification is a convenient
%     method when multiple values need to be input. The user can only change
%     some properties, while leaving others being in their default values.
%     With this function, the function composer can design the function
%     interfaces in the following form: func(main_parameters, varargin), and
%     use parseprops(P0, varargin{:}) to parse the properties specified in the
%     extra arguments.
%         
% $ History $
%   - Created by Dahua Lin on Dec 19th, 2005
%   - Modified by Dahua Lin, on Aug 27, 2006
%       - enhance the compatibility with older version of Matlab
%

if nargin == 1
    P = P0;
elseif nargin == 2
    if isempty(varargin{1})
        P = P0;
    elseif isstruct(varargin{1})
        P = parse_props_from_struct(P0, varargin{1});
    else
        error('sltoolbox:invalidarg', ...
            'Invalid input arguments for properties specification');
    end
else
    n = length(varargin);
    nitems = n / 2;
    if nitems ~= floor(nitems)  % not an integer
        error('sltoolbox:invalidarg', ...
            'Invalid input arguments for properties specification');
    end
    args = reshape(varargin, [2, nitems]);
    P = parse_props_from_arglist(P0, args);
end


%%------------ parse functions  -----------------

function P = parse_props_from_struct(P0, S)

s_names = fieldnames(S);
check_names(P0, s_names);

n = length(s_names);
P = P0;
for i = 1 : n
    curname = s_names{i};
    P.(curname) = S.(curname);
end


function P = parse_props_from_arglist(P0, args)

a_names = args(1, :);
check_names(P0, a_names);

n = length(a_names);
P = P0;
for i = 1 : n
    curname = a_names{i};
    P.(curname) = args{2, i};
end



function check_names(P0, names)

n = length(names);
for i = 1 : n
    if ~ischar(names{i})
        error('Encounter a non-char property name');
    end
end

for i = 1 : n
    if ~isfield(P0, names{i})
        error('The property name %s is invalid', names{i});
    end
end






    
    
    


