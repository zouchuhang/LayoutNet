function dfupdateallplots(dods,dofit,force)
%DFUPDATEALLPLOTS Call update methods for all data sets and fits

%   $Revision: 1.1.6.5 $  $Date: 2004/01/24 09:35:59 $
%   Copyright 2003-2004 The MathWorks, Inc.

le = lasterr;
msg = '';
if nargin<3
   force = false;
end

% Supply defaults if missing or if called as a listener (args not logical)
if nargin<1 || ~islogical(dods)
   dods = true;
end
if nargin<2 || ~islogical(dofit)
   dofit = true;
end

% Get the array of data sets and update each one
if dods
   dsdb = getdsdb;
   ds = down(dsdb);
   while(~isempty(ds))
      try
         if (force)
             clearplot(ds);
         end    
         updateplot(ds);
      catch
         msg = appendmsg(msg,ds.name,lasterr);
      end
      ds = right(ds);
   end
end

% Get the array of fits and update each one
if dofit
   fitdb = getfitdb;
   ft = down(fitdb);
   while(~isempty(ft))
      try
         updateplot(ft);
      catch
         msg = appendmsg(msg,ft.name,lasterr);
      end
      ft = right(ft);
   end
end

lasterr(le);
if ~isempty(msg)
   errordlg(msg,'Error Updating Plot','modal');
end


%--------------------------------------------
function msg = appendmsg(msg,objname,newmsg)
%APPENDMSG Append a new section to an existing set of error messages

if isempty(msg)
   msg = sprintf('Error plotting %s:\n%s',objname,newmsg);
elseif ~isempty(newmsg)
   msg = sprintf('%s\n\nError plotting %s:\n%s',msg,objname,newmsg);
end
