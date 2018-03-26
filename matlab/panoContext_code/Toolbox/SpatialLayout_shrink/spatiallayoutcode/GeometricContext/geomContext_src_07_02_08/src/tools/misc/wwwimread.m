 function [I] = wwwimread(url)

 %create a temporary file name with legal characters
 if 0
 tmpdir = '/tmp/';
 tmpfile = url;
 tmpfile = strrep(tmpfile,'//','');
 tmpfile = strrep(tmpfile,'/','-');
 tmpfile = strrep(tmpfile,'.','-');
 tmpfile = strrep(tmpfile,'http:',tmpdir);%run wget
 command = sprintf('wget %s -O %s',url,tmpfile);
 [s,w] = unix(command);
 end

 url
 s = urlread(url);
 size(s)
 fid = fopen(['./tmp.' url(end-2:end)], 'w');
 fwrite(fid, s);
 fclose(fid);
 
 I = imread(['./tmp.' url(end-2:end)]);
 delete(['./tmp.' url(end-2:end)]);

 if 0
 if (s == 0)
   %read in the file then delete the tmpfile
   I = imread(tmpfile);
   command = sprintf('rm -f %s',tmpfile);
   [s,w] = unix(command);
 else
   disp('error!');
   I = [];
 end;

end