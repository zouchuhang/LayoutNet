function dfhelpviewer(topic, errorname)
% DFHELPVIEWER  is a helper file for the Distribution Fitting Toolbox 
% DFHELPVIEWER Displays help for Distriubtion Fitting TOPIC. If the map file 
% cannot be found, an error is displayed using ERRORNAME

%   Copyright 2003-2004 The MathWorks, Inc. 
%   $Revision: 1.1.6.2 $

import java.io.*;

error = false;
mapfilename = [docroot '/toolbox/stats/stats.map'];
f = File(mapfilename);
if f.exists
    try
        helpview(mapfilename, topic);
    catch
        error = true;
    end
else
    error = true;
end
if error
	message = sprintf('Unable to display help for %s\n', ...
							errorname);
	errordlg(message);
end
