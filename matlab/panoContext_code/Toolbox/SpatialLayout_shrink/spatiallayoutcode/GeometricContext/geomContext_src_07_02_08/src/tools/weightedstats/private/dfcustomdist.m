function dfcustomdist(ignore1,ignore2,action)
%DFCUSTOM Callbacks for menu items related to custom distributions

%   $Revision: 1.1.6.3 $  $Date: 2004/02/01 22:10:36 $
%   Copyright 2003-2004 The MathWorks, Inc.

fnpath = which('dfittooldists.m');
dft = com.mathworks.toolbox.stats.DistributionFitting.getDistributionFitting;

switch(action)
  % --------------------------------
  case 'clear'        % clear custom definitions, revert to standard ones
   % Ask for confirmation
   ok = questdlg('Clear custom distributions and revert to standard ones?',...
                 'Clear Custom Distributions',...
                 'Yes','No','Yes');
   if ~isequal(ok,'Yes')
      return;
   end
   dfgetset('alldistributions','');       % clear all distributions
   dists = dfgetdistributions('',false);  % get built-in list
   dfsetdistributions(dft,dists);         % save them as current
   showresults({dists.name},[]);

  % --------------------------------
  case 'define'       % define a file of probability distribution specs
   % Determine if such a file already exists
   if isempty(fnpath)
      % None found, so start editing a new one with default contents
      fname = fullfile(matlabroot,'toolbox','stats','private',...
                       'dftoolinittemplate.m');
      txt = textread(fname,'%s','whitespace','','bufsize',1e6);
      com.mathworks.mlservices.MLEditorServices.newDocument(txt{1});
   else
      % Edit the existing file
      edit(fnpath)
   end
   
   % Display a helpful message about what's going on
   msg = sprintf(['Define your custom distributions by editing this file and'...
                  '\nsaving it on your path with the name dfittooldists.m.'...
                  '\n\nThen use File -> Custom Distributions -> Import '...
                  '\nto import your distributions.']);
   msgbox(msg,'Define Custom Distributions','none','modal');

  % --------------------------------
  case 'import'       % import a file of probability distribution specs 
   % Remember current distributions
   olds = dfgetset('alldistributions');
   
   % Locate the file of new distribution settings
   if isempty(fnpath)
      fnpath = '*.m';
   end
   [fn,pn] = uigetfile(fnpath,'Select file of distributions to import');
   if isequal(fn,0)
      return
   end
   [dirpath,fname,fext] = fileparts(fn);
   if ~isequal(fext,'.m')
      errordlg(sprintf(['MATLAB .m file required.\n' ...
                        'Can''t import distributions from the file %s.'],...
                       [fname fext]),...
               'Bad Selection');
      return
   end

   % Go to that file's directory and try to run it
   olddir = pwd;
   dists = olds;
   try
      cd(pn);
      [dists,errid,errmsg,newrows] = dfgetuserdists(dists,fname);
   catch
      errmsg = lasterr;
      newrows = [];
   end
   cd(olddir);
   
   % Revert to previous distribution list if anything bad happened
   if ~isempty(errmsg)
      dists = olds;
      errordlg(sprintf('Error trying to import custom distributions:\n%s',...
                       errmsg),...
               'Import Custom Distributions','modal');
      newrows = [];
   end

   % Sort by name
   lowernames = lower(strvcat(dists.name));
   [ignore, ind] = sortrows(lowernames);
   dists = dists(ind);
   newrows = find(ismember(ind,newrows));

   if isempty(errmsg)
      showresults({dists.name},newrows);
   end
   dfsetdistributions(dft,dists);
end

% ---------------------------------
function showresults(liststring,asterisk);
%SHOWRESULTS Stripped-down version of listdlg, just to show a list

promptstring = 'New parametric distribution list:';

if nargin>=2
   for j=1:length(asterisk)
      liststring{asterisk(j)} = sprintf('%s *',liststring{asterisk(j)});
   end
   footnote = ~isempty(asterisk);
else
   footnote = false;
end

ex = get(0,'defaultuicontrolfontsize')*1.7;  % extent height per line
fp = get(0,'defaultfigureposition');
fus = 8;       % frame/uicontrol spacing
ffs = 8;       % frame/figure spacing
uh = 22;       % uicontrol button height
listsize = [160 300];
if footnote
   footnoteheight = 2*ex;
else
   footnoteheight = 0;
end

w = 2*(fus+ffs)+listsize(1);
h = 2*ffs+6*fus+ex+listsize(2)+uh + footnoteheight;
fp = [fp(1) fp(2)+fp(4)-h w h];  % keep upper left corner fixed

figcol = get(0,'defaultUicontrolBackgroundColor');
fig_props = { ...
    'name'                   'Imported Distributions'  ...
    'color'                  figcol ...
    'resize'                 'off' ...
    'numbertitle'            'off' ...
    'menubar'                'none' ...
    'windowstyle'            'modal' ...
    'visible'                'off' ...
    'integerhandle'          'off'    ...
    'handlevisibility'       'callback' ...
    'position'               fp   ...
    'closerequestfcn'        'delete(gcbf)' ...
    'Dock'                   'off' ...
            };
fig = figure(fig_props{:});

posn = [ffs+fus     fp(4)-(ffs+fus+ex) ...
        listsize(1) ex];

uicontrol('style','text','string',promptstring,...
          'horizontalalignment','left','position',posn);

btn_wid = (fp(3)-2*(ffs+fus)-fus)/2;
liststring=cellstr(liststring);
listbox = uicontrol('style','listbox',...
                    'position',[ffs+fus ffs+uh+4*fus+footnoteheight listsize],...
                    'string',liststring,...
                    'backgroundcolor',figcol,...
                    'max',2,...
                    'tag','listbox',...
                    'value',[],...
                    'callback', 'set(gcbo,''value'',0)');

%frameh = uicontrol('style','frame',...
%                   'position',[ffs+fus-1 ffs+fus-1 btn_wid+2 uh+2],...
%                   'backgroundcolor','k');
if footnote
   uicontrol('style','text','string','* Imported or changed',...
             'horizontalalignment','left',...
             'position',[ffs+fus, ffs+fus+uh+footnoteheight/4, listsize(1), footnoteheight]);
end


ok_btn = uicontrol('style','pushbutton',...
                   'string','OK',...
                   'position',[ffs+fus+listsize(1)/2-btn_wid/2 ffs+fus btn_wid uh],...
                   'callback','delete(gcbf)');

% make sure we are on screen
placetitlebar(fig)
set(fig, 'visible','on');
