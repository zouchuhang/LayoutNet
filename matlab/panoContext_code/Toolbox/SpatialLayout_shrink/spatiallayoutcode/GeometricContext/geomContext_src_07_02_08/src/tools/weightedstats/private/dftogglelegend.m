function togglelegend(dffig)
%TOGGLELEGEND Toggle curve fit plot legend on or off

%   $Revision: 1.1.6.2 $  $Date: 2004/01/24 09:35:55 $
%   Copyright 2003-2004 The MathWorks, Inc.

% Get new state -- note uimenu state reflects old state, and
% uitoggletool state reflects new state
h = gcbo;
if isequal(get(h,'Type'),'uimenu')
   onoff = on2off(get(h,'Checked'));
else
   onoff = get(h,'State');
end
dfgetset('showlegend',onoff);

% Change menu state
h = findall(dffig,'Type','uimenu','Tag','showlegend');
if ~isempty(h), set(h,'Checked',onoff); end

% Change button state
h = findall(dffig,'Type','uitoggletool','Tag','showlegend');
if ~isempty(h), set(h,'State',onoff); end

% Forget previous legend location
if isequal(onoff,'off')
   dfgetset('legendpos',[]);
   dfgetset('rlegendpos',[]);
end

% Change legend state
dfupdatelegend(dffig);
