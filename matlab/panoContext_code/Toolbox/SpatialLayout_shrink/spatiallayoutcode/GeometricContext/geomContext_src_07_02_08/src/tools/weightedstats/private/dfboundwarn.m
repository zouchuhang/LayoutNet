function dfboundwarn(h)
%DFBOUNDWARN Helper function for dfittool to warn about conf bounds

%   Copyright 2003-2004 The MathWorks, Inc. 
%   $Revision: 1.1.6.3 $  $Date: 2004/03/02 21:49:23 $

if isa(h,'stats.dfdata')
   % Warn if the function type does not allow data set bounds
   ftype = dfgetset('ftype');
   if isequal(ftype,'pdf')
      pname = 'Density/histogram';
   elseif isequal(ftype,'probplot')
      pname = 'Probability';
   elseif isequal(ftype,'icdf')
      pname = 'Quantile (inverse cdf)';
   else
      return;
   end
   wmsg = sprintf(['%s plots do not include confidence bounds for data sets.\n' ...
                   'Bounds for the data set "%s" may be visible if ' ...
                   'you change the display type.'], ...
                  pname,h.name);
   warndlg(wmsg,'Distribution Fitting Warning','modal');

else % isa(h,'stats.dffit')
   % Warn if the function type does not allow fit bounds
   if ~dfgetset('dobounds')
      ftype = dfgetset('ftype');
      if isequal(ftype,'pdf')
         pname = 'Density';
      elseif isequal(ftype,'probplot')
         pname = 'Probability';
      else
         return;
      end
      wmsg = sprintf(['%s plots do not include confidence bounds.\n' ...
                      'Bounds for the fit "%s" may be visible if ' ...
                      'you change the display type.'], ...
                     pname,h.name);
      warndlg(wmsg,'Distribution Fitting Warning','modal');
   end
end