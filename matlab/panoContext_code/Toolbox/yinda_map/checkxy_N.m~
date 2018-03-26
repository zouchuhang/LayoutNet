function checkxy(x, y, function_name, x_var_name, y_var_name, x_pos, y_pos)
%CHECKXY Check validity of map x and y vectors
%
%   CHECKXY(X, Y, FUNCTION_NAME, X_VAR_NAME, Y_VAR_NAME, X_POS, X_POS)
%   ensures that X and Y are real vectors of matching size and equal NaN
%   locations.

% Copyright 2006-2011 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2011/05/17 02:12:40 $

% Input arguments are not checked for validity.

if ~isempty(x) || ~isempty(y)
    % Check numeric, 2d, and real.
    validateattributes(x, ...
        {'numeric'}, {'real','2d','vector'}, function_name, x_var_name, x_pos);
    validateattributes(y, ...
        {'numeric'}, {'real','2d','vector'}, function_name, y_var_name, y_pos);
    
    if ~isequal(isnan(x), isnan(y))
        error(sprintf('map:%s:inconsistentXY', function_name), ...
            'Function %s expected its %s and %s input arguments, %s and %s, to match in size or NaN locations.', ...
            upper(function_name), num2ordinal(x_pos), num2ordinal(y_pos), ...
            x_var_name, y_var_name)
    end
end
