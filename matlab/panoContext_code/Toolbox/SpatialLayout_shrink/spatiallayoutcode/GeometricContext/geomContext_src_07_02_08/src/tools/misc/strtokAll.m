function tok = strtokAll(str, del)
% tok = strtokAll(str, del)
% Returns a cell array of tokens from str deliminated by del
% example: tok = strtokAll('My cat ate you', ' ') ==>
%               {'My','cat','ate','you'}

tok = {};

rem = str;
while ~isempty(rem)
    [tok{end+1}, rem] = strtok(rem, del);
end


