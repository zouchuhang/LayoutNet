function dffig2m(dffig,outfilename)
%DFFIG2M Turn figure into an M file that can produce the figure

%   $Revision: 1.1.6.12 $  $Date: 2004/03/22 23:55:33 $
%   Copyright 2003-2004 The MathWorks, Inc.

dsdb = getdsdb;
fitdb = getfitdb;
if isempty(down(dsdb)) && isempty(down(fitdb))
   emsg = 'Cannot save M file when no datasets or fits exist.';
   errordlg(emsg,'Error Saving M File','modal');
   return
end

if nargin<1
   dffig = dfgetset('dffig'); 
end

if nargin<2
   % Get file name to use, remember the directory name
   olddir = dfgetset('dirname');
   filespec = [olddir '*.m'];
   [outfilename,pn] = uiputfile(filespec,'Save M File');
   if isequal(outfilename,0) || isequal(pn,0)
      return
   end
   if ~ismember('.',outfilename)
      outfilename = [outfilename '.m'];
   end
   dfgetset('dirname',pn);
   outfilename = sprintf('%s%s',pn,outfilename);
end

% Get M file name with .m suffix, and get corresponding function name
if length(outfilename)<2 || ~isequal(outfilename(end-1:end),'.m')
   outfilename = sprintf('%s.m',outfilename);
end
fcnname = outfilename(1:end-2);
k = max(find(fcnname(1:end-1)=='\'));
if ~isempty(k)
   fcnname = fcnname(k+1:end);
end
k = max(find(fcnname(1:end-1)=='/'));
if ~isempty(k)
   fcnname = fcnname(k+1:end);
end
   

% Set up some variables for later
allprop = {'Color' 'Marker' 'LineStyle' 'LineWidth' 'MarkerSize'};
showlegend = isequal(dfgetset('showlegend'),'on');
ftype = dfgetset('ftype');
alpha = 1 - dfgetset('conflev');

% Create arrays to receive code text
blkc = cell(0,1);    % block of comment lines
blks = cell(0,1);    % block of setup lines
blkd = cell(0,1);    % block of data-related lines
blkf = cell(0,1);    % block of fit-related lines
blke = cell(0,1);    % block of lines at end

% Write introduction to dataset section, including figure
% preparation code
blks{end+1} = '% Set up figure to receive datasets and fits';
blks{end+1} = 'f_ = clf;';
blks{end+1} = 'figure(f_);';
if showlegend
   blks{end+1} = 'legh_ = []; legt_ = {};   % handles and text for legend';
end

% Process each dataset
exprlist = {};    % names and expressions of the data, censoring, frequency
arglist = {};     % variable names to use for each expression
ds = down(dsdb);
numds = 0;
while(~isempty(ds))
   numds = numds + 1;
   [blkc,blkd,exprlist,arglist,showbounds,onplot] = ...
                 writedset(blkc,blkd,ds,exprlist,arglist,allprop,alpha);

   if onplot && showlegend
      blkd{end+1} = 'legh_(end+1) = h_;';
      blkd{end+1} = sprintf('legt_{end+1} = ''%s'';',quotedtext(ds.name));
      if showbounds
         blkd{end+1} = 'legh_(end+1) = hb_;';
         blkd{end+1} = sprintf('legt_{end+1} = ''%g%% confidence bounds'';',...
                               100*(1-alpha));
      end
   end
   ds = right(ds);
end

% Set up for plotting fits
anycontinuous = false;
anydiscrete = false;
ft = down(fitdb);
while(~isempty(ft))
   if ft.iscontinuous
      anycontinuous = true;
   else
      anydiscrete = true;
   end
   ft = right(ft);
end

% Create a suitable X vector, may depend on whether it's discrete
if ~isequal(ftype,'pdf') || ~anydiscrete
   blkf{end+1} = sprintf('x_ = linspace(xlim_(1),xlim_(2),100);');
elseif ~anycontinuous
   blkf{end+1} = 'incr_ = max(1,floor((xlim_(2)-xlim_(1))/100));';
   blkf{end+1} = 'x_ = floor(xlim_(1)):incr_:ceil(xlim_(2));';
else
   blkf{end+1} = sprintf('xc_ = linspace(xlim_(1),xlim_(2),100);');
   blkf{end+1} = 'incr_ = max(1,floor((xlim_(2)-xlim_(1))/100));'
   blkf{end+1} = 'xd_ = floor(xlim_(1)):incr_:ceil(xlim_(2));';
end

% Process each fit
numfit = 0;
ft = down(fitdb);
anySmoothFits = false;
while(~isempty(ft))
   numfit = numfit+1;
   fitname = ft.name;

   % Create code to re-create this fit
   blkf{end+1} = sprintf('\n%% --- Create fit "%s"',fitname);

   % Call subfunction to generate code for each type
   if isequal(getfittype(ft),'param')
      [blkf,showbounds,onplot] = writepfit(blkf,ft,alpha,allprop,...
                                    anycontinuous,anydiscrete,exprlist,arglist);
   else
      anySmoothFits = true;
      [blkf,onplot] = writenpfit(blkf,ft,alpha,allprop,...
                                 anycontinuous,anydiscrete,exprlist,arglist);
      showbounds = false;
   end

   % Add legend if requested
   if onplot && showlegend
      blkf{end+1} = 'legh_(end+1) = h_;';
      blkf{end+1} = sprintf('legt_{end+1} = ''%s'';',quotedtext(ft.name));
      if showbounds
         blkf{end+1} = 'legh_(end+1) = hb_;';
         blkf{end+1} = sprintf('legt_{end+1} = ''%g%% confidence bounds'';',...
                               100*(1-alpha));
      end
   end
   ft = right(ft);
end

% In setup section, create empty axes and set some properties
if ~isequal(ftype,'probplot')
   blks{end+1} = 'ax_ = newplot;';
else
   dtype = dfgetset('dtype');
   if ischar(dtype)
      blks{end+1} = sprintf('probplot(''%s'');', dtype);
   else
      blks{end+1} = sprintf(...
          'dist_ = dfswitchyard(''dfgetdistributions'',''%s'');',...
          dtype.distspec.code);
      blks{end+1} = sprintf('probplot({dist_,%s});',...
                            cell2text(num2cell(dtype.params)));
   end
   blks{end+1} = 'ax_ = gca;';
   blks{end+1} = 'title(ax_,'''');';
end

blks{end+1} = 'set(ax_,''Box'',''on'');';
if isequal(dfgetset('showgrid'),'on')
   blks{end+1} = 'grid(ax_,''on'');';
end
blks{end+1} = 'hold on;';

% At end of data set section, set x axis limits
blkd{end+1} = sprintf('\n%% Nudge axis limits beyond data limits');
blkd{end+1} = 'xlim_ = get(ax_,''XLim'');';
blkd{end+1} = 'if all(isfinite(xlim_))';
blkd{end+1} = '   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);';
blkd{end+1} = '   set(ax_,''XLim'',xlim_)';
blkd{end+1} = 'end';


% Finish up
blke{end+1} = 'hold off;';
if showlegend
   axold = get(dffig,'CurrentAxes');
   legh = legend('-find',axold);
   if isempty(legh) || ~ishandle(legh)
      oldpos = [.8 .8 .1 .1];  % close to TR position
   else
      oldpos = get(legh,'Position');
   end
   if oldpos(1)<.4
      newpos = 'NorthWest';
   else
      newpos = 'NorthEast';
   end

   blke{end+1} = sprintf('legend(ax_,legh_, legt_, ''Location'',''%s'');',newpos);
end

% Write code into m file
if length(arglist)==0
   argtext = '';
else
   argtext = sprintf('%s,',arglist{:});
   argtext = sprintf('(%s)',argtext(1:end-1));
end
[fid,msg] = fopen(outfilename,'w');
if fid==-1
   emsg = sprintf('Error trying to write to %s:\n%s',outfilename,msg);
   errordlg(emsg,'Error Saving M File','modal');
   return
end
fprintf(fid,'function %s%s\n',fcnname,argtext);
fprintf(fid,'%%%s    Create plot of datasets and fits\n',upper(fcnname));
fprintf(fid,'%%   %s%s\n',upper(fcnname),upper(argtext));
fprintf(fid,'%%   Creates a plot, similar to the plot in the main distribution fitting\n');
fprintf(fid,'%%   window, using the data that you provide as input.  You can\n');
fprintf(fid,'%%   apply this function to the same data you used with dfittool\n');
fprintf(fid,'%%   or with different data.  You may want to edit the function to\n');
fprintf(fid,'%%   customize the code and this help message.\n');
fprintf(fid,'%%\n');
fprintf(fid,'%%   Number of datasets:  %d\n',numds);
fprintf(fid,'%%   Number of fits:  %d\n',numfit);
fprintf(fid,'\n');
fprintf(fid,'%% This function was automatically generated on %s\n',...
            datestr(now));
for j=1:length(blkc)
   fprintf(fid,'%s\n',xlate(blkc{j}));
end
fprintf(fid,'\n');
for j=1:length(blks)
   fprintf(fid,'%s\n',xlate(blks{j}));
end
fprintf(fid,'\n');
for j=1:length(blkd)
   fprintf(fid,'%s\n',xlate(blkd{j}));
end
fprintf(fid,'\n');
for j=1:length(blkf)
   fprintf(fid,'%s\n',xlate(blkf{j}));
end
fprintf(fid,'\n');
for j=1:length(blke)
   fprintf(fid,'%s\n',xlate(blke{j}));
end

% Create sub function to be used to support a functionline fit on a probability plot
if anySmoothFits && isequal(ftype,'probplot')
   fprintf(fid,'\n\n%% -----------------------------------------------\n');
   fprintf(fid,'function f=cdfnp(x,y,cens,freq,support,kernel,width)\n');
   fprintf(fid,'%%CDFNP Compute cdf for non-parametric fit, used in probability plot\n\n');
   fprintf(fid,'f = ksdensity(y,x,''cens'',cens,''weight'',freq,''function'',''cdf'',...\n');
   fprintf(fid,'                  ''support'',support,''kernel'',kernel,''width'',width);\n');
end

fclose(fid);

% ------------------- double up quotes in text string
function a = quotedtext(b)
if ischar(b)
   a = strrep(b,'''','''''');
else
   a = sprintf('%.13g',b);
end

% ------------------- create text to re-create cell or numeric array
function a = cell2text(b)

if ~iscell(b)
   if ischar(b)
      a = sprintf('''%s''',quotedtext(b));
   elseif length(b)==1
      a = sprintf('%.13g',b);
   else
      numtext = num2str(b,'%.13g ');
      if size(numtext,1)>1
         numtext = [numtext repmat(';',size(numtext,1),1)]';
         numtext = numtext(:)';
         numtext = numtext(1:end-1);
      end
      a = sprintf('[%s]',numtext);
   end
   return
end

if length(b)>0
   bj = b{1};
   if ischar(bj)
      a = sprintf('''%s''',quotedtext(bj));
   else
      a = sprintf('%.13g',bj);
   end
   for j=2:length(b)
      bj = b{j};
      if ischar(bj)
         a = sprintf('%s, ''%s''',a,quotedtext(bj));
      else
         a = sprintf('%s, %.13g',a,bj);
      end
   end
else
   a = '';
end
a = sprintf('[%s]',a);


% ----------------- add censoring and frequency args to code block
function blk = addcensfreq(blk,censname,freqname)

if ~isempty(censname) && ~isequal(censname,'[]')
   blk{end+1} = sprintf('               ,''cens'',%s...',censname);
end
if ~isempty(freqname) && ~isequal(freqname,'[]')
   blk{end+1} = sprintf('               ,''freq'',%s...',freqname);
end


% ---------------- write code for parametric fit
function [blkf,showbounds,onplot] = ...
    writepfit(blkf,ft,alpha,allprop,anycontinuous,anydiscrete,exprlist,arglist)

ds = ft.ds;
yname = expression2name(ds.yname,exprlist,arglist);
dist = ft.distspec;
ftype = ft.ftype;
showbounds = false;
onplot = true;

blkf{end+1} = sprintf('\n%% Fit this distribution to get parameter values');
[censname,freqname] = getcensfreqname(ds,exprlist,arglist);
shortform = isempty(censname) & isempty(freqname);

% Exclude data if necessary
if ~isempty(ft.exclusionrule)
   [blkf,yname,censname,freqname] = applyexclusion(blkf,ft.exclusionrule,...
                                                   yname,censname,freqname);
end

if isempty(censname)
   censname = '[]';
end
if isempty(freqname)
   freqname = '[]';
end

% Helpful note about using old results instead of fitting new data
if isequal(getfittype(ft),'param')
    blkf{end+1} = sprintf('%% To use parameter estimates from the original fit:');
    blkf{end+1} = sprintf('%%     p_ = %s;', cell2text(num2cell(ft.params)));
end

nparams = length(dist.pnames);

if shortform
   arglist = sprintf('%s, %g',yname,alpha);
else
   arglist = sprintf('%s, %g, %s, %s',yname,alpha,censname,freqname);
end

fname = func2str(dist.fitfunc);
onpath = exist(fname);
if onpath
   rhs = sprintf('%s(%s);',fname,arglist);
else
   rhs = sprintf('mle(''%s'',%s);  %% Fit %s distribution',...
                 dist.code,arglist,dist.name);
end
if dist.paramvec
   blkf{end+1} = sprintf('p_ = %s;',rhs);
else
   blkf{end+1} = sprintf('pargs_ = cell(1,%d);',nparams);
   blkf{end+1} = sprintf('[pargs_{:}] = %s',rhs);
   blkf{end+1} = 'p_ = [pargs_{:}];';
end

pargs = 'p_(1)';
if nparams>1
   pargs = [pargs, sprintf(', p_(%d)',2:nparams)];
end

% Get covariance matrix if we need confidence bounds
if ft.showbounds && ismember(ftype,{'cdf' 'survivor' 'cumhazard' 'icdf'})
   showbounds = true;
else
   showbounds = false;
end

% Sometimes we need the structure that describes the distribution
if ~onpath && (showbounds || isequal(ftype,'probplot'))
   blkf{end+1} = sprintf(...
          '\n%% Get a description of the %s distribution',...
          dist.name);
   blkf{end+1} = sprintf(...
          'dist_ = dfswitchyard(''dfgetdistributions'',''%s'');\n',...
          dist.code);
end

if showbounds
   if onpath
      blkf{end+1} = sprintf('[NlogL_,pcov_] = %s(p_,%s,%s,%s);',...
                            func2str(dist.likefunc),yname, censname, freqname);
   else
      blkf{end+1} = sprintf(...
          '[NlogL_,pcov_] = feval(dist_.likefunc,p_,%s,%s,%s);',...
          yname, censname, freqname);
   end
end

% Plot the fit and bounds if the original figure had them plotted
if isempty(ft.line) || ~ishandle(ft.line)
   blkf{end+1} = '% This fit does not appear on the plot';
   onplot = false;
   return;
end

propvals = get(ft.line,allprop);
[c,m,l,w,s] = deal(propvals{:});

switch(ftype)
 case {'pdf'}
   if anycontinuous && anydiscrete
      if ft.iscontinuous
         blkf{end+1} = 'x_ = xc_;';
      else
         blkf{end+1} = 'x_ = xd_;';
      end            
   end
   if onpath
      blkf{end+1} = sprintf('y_ = %s(x_,%s);',func2str(dist.pdffunc),pargs);
   else
      blkf{end+1} = sprintf('y_ = pdf(''%s'',x_,%s);',dist.code,pargs);
   end
   
   blkf{end+1} = sprintf('h_ = plot(x_,y_,''Color'',[%g %g %g],...',...
                         c(1),c(2),c(3));
   blkf{end+1} = sprintf('          ''LineStyle'',''%s'', ''LineWidth'',%d,...',l,w);
   blkf{end+1} = sprintf('          ''Marker'',''%s'', ''MarkerSize'',%d);',m,s);
  
 case {'cdf' 'survivor' 'cumhazard' 'icdf'}
   if isequal(ftype,'icdf')
      if onpath
         prefix = sprintf('%s(',func2str(dist.invfunc));
      else
         prefix = sprintf('icdf(''%s'',',dist.code);
      end
   else
      if onpath
         prefix = sprintf('%s(',func2str(dist.cdffunc));
      else
         prefix = sprintf('cdf(''%s'',',dist.code);
      end
   end

   if showbounds
      blkf{end+1} = sprintf('[y_,yL_,yU_] = %sx_,%s,pcov_,%g); %% cdf and bounds',...
                            prefix,pargs,alpha);
   else
      blkf{end+1} = sprintf('y_ = %sx_,%s); %% compute cdf',...
                            prefix,pargs);
   end

   if isequal(ftype,'survivor')
      blkf{end+1} = 'y_ = 1 - y_; % convert to survivor function';
      if showbounds
         blkf{end+1} = 'tmp_ = yL_;';
         blkf{end+1} = 'yL_ = 1 - yU_;';
         blkf{end+1} = 'yU_ = 1 - tmp_;';
      end
   elseif isequal(ftype,'cumhazard')
      blkf{end+1} = 't_ = (y_ < 1); % only where the hazard is finite';
      blkf{end+1} = 'x_ = x_(t_);';
      blkf{end+1} = 'y_ = -log(1 - y_(t_));';
      if showbounds
         blkf{end+1} = 'if ~isempty(yL_)';
         blkf{end+1} = '   tmp_ = yL_;';
         blkf{end+1} = '   yL_ = -log(1 - yU_(t_));';
         blkf{end+1} = '   yU_ = -log(1 - tmp_(t_));';
         blkf{end+1} = 'end';
      end
   end
      
   blkf{end+1} = sprintf('h_ = plot(x_,y_,''Color'',[%g %g %g],...',...
                         c(1),c(2),c(3));
   blkf{end+1} = sprintf('          ''LineStyle'',''%s'', ''LineWidth'',%d,...',l,w);
   blkf{end+1} = sprintf('          ''Marker'',''%s'', ''MarkerSize'',%d);',m,s);

   if showbounds
      blkf{end+1} = 'if ~isempty(yL_)';
      blkf{end+1} = sprintf('   hb_ = plot([x_(:); NaN; x_(:)], [yL_(:); NaN; yU_(:)],''Color'',[%g %g %g],...',...
                            c(1),c(2),c(3));
      blkf{end+1} = '             ''LineStyle'','':'', ''LineWidth'',1,...';
      blkf{end+1} = '             ''Marker'',''none'');';
      blkf{end+1} = 'end';
   end

 case 'probplot'
   if onpath
      stmt = sprintf('h_ = probplot(ax_,@%s,p_);', ...
                     func2str(dist.cdffunc));
   else
      stmt = 'h_ = probplot(ax_,dist_.cdffunc,p_);';
   end
   blkf{end+1} = stmt;
   blkf{end+1} = sprintf('set(h_,''Color'',[%g %g %g],''LineStyle'',''%s'', ''LineWidth'',%d);', ...
         c(1),c(2),c(3),l,w);
end


% ---------------- write code for nonparametric fit
function [blkf,onplot] = ...
    writenpfit(blkf,ft,alpha,allprop,anycontinuous,anydiscrete,exprlist,arglist)

ds = ft.ds;
yname = expression2name(ds.yname,exprlist,arglist);
ftype = ft.ftype;

[censname,freqname] = getcensfreqname(ds,exprlist,arglist);
shortform = isempty(censname) & isempty(freqname);
   
% Exclude data if necessary
if ~isempty(ft.exclusionrule)
   [blkf,yname,censname,freqname] = applyexclusion(blkf,ft.exclusionrule,...
                                                   yname,censname,freqname);
end

if isempty(censname)
   censname = '[]';
end
if isempty(freqname)
   freqname = '[]';
end

kernel = sprintf('''%s''',ft.kernel);
if ft.bandwidthradio == 0
   width = '[]';
else
   width = ft.bandwidthtext;
end
if ischar(ft.support)
   spt = sprintf('''%s''',ft.support);
else
   spt = sprintf('[%g, %g]',ft.support);
end

% Plot the fit and bounds if the original figure had them plotted
if isempty(ft.line) || ~ishandle(ft.line)
   blkf{end+1} = '% This fit does not appear on the plot';
   onplot = false;
   return;
end
onplot = true;

propvals = get(ft.line,allprop);
[c,m,l,w,s] = deal(propvals{:});

switch(ftype)
 case {'pdf' 'icdf' 'cdf' 'survivor' 'cumhazard'}
   if isequal(ftype,'pdf') && anycontinuous && anydiscrete
      blkf{end+1} = 'x_ = xc_;';
   end

   blkf{end+1} = sprintf('y_ = ksdensity(%s,x_,''kernel'',%s,...',...
                         yname,kernel);
   if ~shortform
      blkf{end+1} = sprintf('               ''cens'',%s,''weight'',%s,...',...
                         censname,freqname);
   end

   if ~isequal(ft,'unbounded')
      blkf{end+1} = sprintf('               ''support'',%s,...',spt);
   end
   if ~isequal(width,'[]')
      blkf{end+1} = sprintf('               ''width'',%s,...',width);
   end      
   
   blkf{end+1} = sprintf('               ''function'',''%s'');',ftype);
   blkf{end+1} = sprintf('h_ = plot(x_,y_,''Color'',[%g %g %g],...',...
                         c(1),c(2),c(3));
   blkf{end+1} = sprintf('          ''LineStyle'',''%s'', ''LineWidth'',%d,...',l,w);
   blkf{end+1} = sprintf('          ''Marker'',''%s'', ''MarkerSize'',%d);',m,s);
  

 case 'probplot'
   blkf{end+1} = sprintf('npinfo_ = {%s %s %s %s %s %s};',...
                         yname,censname,freqname,spt,kernel,width);

   blkf{end+1} = 'h_ = probplot(ax_,@cdfnp,npinfo_);';
   blkf{end+1} = sprintf('set(h_,''Color'',[%g %g %g],''LineStyle'',''%s'', ''LineWidth'',%d);', ...
         c(1),c(2),c(3),l,w);
end


% --------------- write code for data set
function [blkc,blkd,exprlist,arglist,showbounds,onplot] = ...
                 writedset(blkc,blkd,ds,exprlist,arglist,allprop,alpha)

dsname = ds.name;
yname = ds.yname;
[censname,freqname] = getcensfreqname(ds);
newnames = {yname censname freqname};
newvars = false(1,3);
ftype = ds.ftype;
showbounds = false;
onplot = true;

% Create comment text associating dataset with variable names
blkc{end+1} = ' ';
blkc{end+1} = sprintf('%% Data from dataset "%s":',dsname);

% Each non-empty variable name becomes a function argument,
% except expressions that are not valid variable names have
% to be replaced by a variable name that we will select here
descrtext = {'Y' 'Censoring' 'Frequency'};
for j=1:3
   exprj = newnames{j};
   if isempty(exprj)
      continue;
   end
   exprnum = strmatch(exprj,exprlist,'exact');
   if isempty(exprnum);
      exprnum = length(exprlist) + 1;
      exprlist{exprnum} = exprj;
      if isvarname(exprj)
         namej = exprj;
      else
         namej = sprintf('arg_%d',exprnum);
      end
      arglist{exprnum} = namej;
      newvars(j) = true;
   else
      namej = arglist{exprnum};
   end
   if isequal(namej,exprj)
      suffix = '';
   else
      suffix = sprintf(' (originally %s)',exprj);
   end

   blkc{end+1} = sprintf('%%    %s = %s%s',descrtext{j},namej,suffix);
   newnames{j} = namej;
end

yname = newnames{1};
censname = newnames{2};
freqname = newnames{3};
havecens = ~isempty(censname);
if ~havecens
   censname = '[]';
end
havefreq = ~isempty(freqname);
if ~havefreq
   freqname = '[]';
end

blkc{end+1} = ' ';
blkc{end+1} = '% Remove missing values';
blkc{end+1} = sprintf('t_ = ~isnan(%s);', yname);
if havecens
   blkc{end+1} = sprintf('if ~isempty(%s), t_ = t_ & ~isnan(%s); end',...
                         censname, censname);
end
if havefreq
   blkc{end+1} = sprintf('if ~isempty(%s), t_ = t_ & ~isnan(%s); end',...
                         freqname, freqname);
end
blkc{end+1} = sprintf('%s = %s(t_);', yname, yname);
if havecens
   blkc{end+1} = sprintf('if ~isempty(%s), %s = %s(t_); end',...
                         censname, censname, censname);
end
if havefreq
   blkc{end+1} = sprintf('if ~isempty(%s), %s = %s(t_); end',...
                         freqname, freqname, freqname);
end


% Create code to plot this dataset into the figure we have created
blkd{end+1} = ' ';
blkd{end+1} = sprintf('%% --- Plot data originally in dataset "%s"',dsname);
for j=1:3
   if newvars(j)
      blkd{end+1} = sprintf('%s = %s(:);',newnames{j},newnames{j});
   end
end
dsline = ds.line;
if isempty(dsline) || ~ishandle(dsline)
   blkd{end+1} = '% This dataset does not appear on the plot';
   onplot = false;
   return;
end

propvals = get(dsline,allprop);
[c,m,l,w,s] = deal(propvals{:});
switch(ftype)
 case 'pdf'
   % Generate code to compute the empirical cdf
   blkd{end+1} = sprintf('[F_,X_] = ecdf(%s,''Function'',''cdf''...', yname);
   if havecens
      blkd{end+1} = sprintf('               ,''cens'',%s...',censname);
   end
   if havefreq
      blkd{end+1} = sprintf('               ,''freq'',%s...',freqname);
   end
   blkd{end+1} = '              );  % compute empirical cdf';

   % Generate code to duplicate the current histogram bin width selection
   bininfo = ds.binDlgInfo;
   if isempty(bininfo)           % use default in case this is empty
      bininfo.rule = 1;
   end
   blkd{end+1} = sprintf('Bin_.rule = %d;', bininfo.rule);

   switch bininfo.rule
    case 3
      blkd{end+1} = sprintf('Bin_.nbins = %d;',bininfo.nbins);

    case 5
      blkd{end+1} = sprintf('Bin_.width = %g;',bininfo.width);
      blkd{end+1} = sprintf('Bin_.placementRule = %d;',bininfo.placementRule);
      if bininfo.placementRule ~= 1
         blkd{end+1} = sprintf('Bin_.anchor = %g;',bininfo.anchor);
      end
   end
   
   blkd{end+1} = sprintf('[C_,E_] = dfswitchyard(''dfhistbins'',%s,%s,%s,Bin_,F_,X_);',...
                         yname,censname,freqname);

   % Generate code to produce the histogram
   blkd{end+1} = '[N_,C_] = ecdfhist(F_,X_,''edges'',E_); % empirical pdf from cdf';
   blkd{end+1} = 'h_ = bar(C_,N_,''hist'');';
   blkd{end+1} = sprintf('set(h_,''FaceColor'',''none'',''EdgeColor'',[%g %g %g],...', ...
                         c(1),c(2),c(3));
   blkd{end+1} = sprintf('       ''LineStyle'',''%s'', ''LineWidth'',%d);', ...
         l,w);
   blkd{end+1} = 'xlabel(''Data'');';
   blkd{end+1} = 'ylabel(''Density'')';
  
 case {'cdf' 'survivor' 'cumhazard'}
   showbounds = ds.showbounds;
   if showbounds
      blkd{end+1} = sprintf('[Y_,X_,yL_,yU_] = ecdf(%s,''Function'',''%s'',''alpha'',%g...',...
                            yname, ftype,alpha);
   else
      blkd{end+1} = sprintf('[Y_,X_] = ecdf(%s,''Function'',''%s''...',...
                            yname, ftype);
   end
   blkd = addcensfreq(blkd,censname,freqname);
   blkd{end+1} = '              );  % compute empirical function';
   blkd{end+1} = 'h_ = stairs(X_,Y_);';
   blkd{end+1} = sprintf('set(h_,''Color'',[%g %g %g],''LineStyle'',''%s'', ''LineWidth'',%d);', ...
         c(1),c(2),c(3),l,w);
   if showbounds
      blkd{end+1} = '[XX1_,YY1_] = stairs(X_,yL_);';
      blkd{end+1} = '[XX2_,YY2_] = stairs(X_,yU_);';
      blkd{end+1} = 'hb_ = plot([XX1_(:); NaN; XX2_(:)], [YY1_(:); NaN; YY2_(:)],...';
      blkd{end+1} = sprintf('   ''Color'',[%g %g %g],''LineStyle'','':'', ''LineWidth'',1);', ...
         c(1),c(2),c(3));
   end
   blkd{end+1} = 'xlabel(''Data'');';
   switch(ftype)
      case 'cdf',       blkd{end+1} = 'ylabel(''Cumulative probability'')';
      case 'survivor',  blkd{end+1} = 'ylabel(''Survivor function'')';
      case 'cumhazard', blkd{end+1} = 'ylabel(''Cumulative hazard'')';
   end

 case 'icdf'
   blkd{end+1} = sprintf('[Y_,X_] = ecdf(%s,''Function'',''cdf''...', yname);
   blkd = addcensfreq(blkd,censname,freqname);
   blkd{end+1} = '              );  % compute empirical cdf';
   blkd{end+1} = 'h_ = stairs(Y_,[X_(2:end);X_(end)]);';
   blkd{end+1} = sprintf('set(h_,''Color'',[%g %g %g],''LineStyle'',''%s'', ''LineWidth'',%d);', ...
         c(1),c(2),c(3),l,w);
   blkd{end+1} = 'xlabel(''Probability'');';
   blkd{end+1} = 'ylabel(''Quantile'')';

 case 'probplot'
   blkd{end+1} = sprintf('h_ = probplot(ax_,%s...', yname);
   blkd = addcensfreq(blkd,censname,freqname);
   blkd{end+1} = '             ,''noref'');  % add to probability plot';
   blkd{end+1} = sprintf('set(h_,''Color'',[%g %g %g],''Marker'',''%s'', ''MarkerSize'',%d);', ...
         c(1),c(2),c(3),m,s);
   blkd{end+1} = 'xlabel(''Data'');';
   blkd{end+1} = 'ylabel(''Probability'')';
end


% -----------------------------
function [blkf,yname,censname,freqname]=applyexclusion(blkf,exclrule,...
                                                       yname,censname,freqname);
%APPLYEXCLUSION Change var names to use indexing to apply exclusion rule

% Create expressions for inclusion rules
if isempty(exclrule.ylow)
   e1 = '';
else
   ylow = str2double(exclrule.ylow);
   if exclrule.ylowlessequal==1
      e1 = sprintf('%s > %g', yname, ylow);
   else
      e1 = sprintf('%s >= %g', yname, ylow);
   end
end
if isempty(exclrule.yhigh)
   e2 = '';
else
   yhigh = str2double(exclrule.yhigh);
   if exclrule.yhighgreaterequal==1
      e2 = sprintf('%s < %g', yname, yhigh);
   else
      e2 = sprintf('%s <= %g', yname, yhigh);
   end
end

% Combine exclusion expressions
if isempty(e1)
   if isempty(e2)
      etxt = '';
   else
      etxt = e2;
   end
else
   if isempty(e2)
      etxt = e1;
   else
      etxt = sprintf('%s & %s',e1,e2);
   end
end

% Create code to generate index vector and reduce all variables
if ~isempty(etxt)
   blkf{end+1} = sprintf('\n%% Create vector for exclusion rule ''%s''',...
                         exclrule.name);
   blkf{end+1} =         '% Vector indexes the points that are included';
   blkf{end+1} = sprintf('excl_ = (%s);\n', etxt);

   yname = sprintf('%s(excl_)',yname);
   if ~isempty(censname)
      censname = sprintf('%s(excl_)',censname);
   end
   if ~isempty(freqname)
      freqname = sprintf('%s(excl_)',freqname);
   end
end

% -----------------------------------------
function [censname,freqname] = getcensfreqname(ds,exprlist,arglist)
%GETCENSFREQNAME Get censoring and freqency names

censname = ds.censname;
freqname = ds.freqname;
if strcmp(censname,'(none)')
   censname = '';
end
if strcmp(freqname,'(none)')
   freqname = '';
end

if nargin>=3
   censname = expression2name(censname,exprlist,arglist);
   freqname = expression2name(freqname,exprlist,arglist);
end
   
  
% -------------------------------------------
function nm = expression2name(expr,exprlist,arglist)
%EXPRESSION2NAME Find out what name we're using in place of this expression

nm = expr;
if ~isempty(expr)
   j = strmatch(expr,exprlist,'exact');
   if isscalar(j)
      nm = arglist{j};
   end
end
   