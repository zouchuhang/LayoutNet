function toks = split_string(line, delim)
ntok = 0;
tmp1 = line;
while length(tmp1)>0
    [tmp2, tmp1] = strtok(tmp1, delim);
    ntok = ntok+1;
end
toks = cell(ntok, 1);
for i = 1:ntok
    [toks{i}, line] = strtok(line, delim);
end