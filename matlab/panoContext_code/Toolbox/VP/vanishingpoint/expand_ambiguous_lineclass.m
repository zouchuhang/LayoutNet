function lines_expand = expand_ambiguous_lineclass(lines)
% when lines belong to more than 1 class,
% hallucinate line for each class.

%%
num_lines = length(lines);

%%
lines_expand(num_lines*2).point1 = [0 0];
lines_expand(num_lines*2).point2 = [0 0];
lines_expand(num_lines*2).lineclass = 0;
count = 0;

%%
for i = 1:num_lines
    if lines(i).lineclass1 == 1
        count = count + 1;
        lines_expand(count).point1 = lines(i).point1;
        lines_expand(count).point2 = lines(i).point2;
        lines_expand(count).lineclass = 1;
    end
    if lines(i).lineclass2 == 1
        count = count + 1;
        lines_expand(count).point1 = lines(i).point1;
        lines_expand(count).point2 = lines(i).point2;
        lines_expand(count).lineclass = 2;
    end
    if lines(i).lineclass3 == 1
        count = count + 1;
        lines_expand(count).point1 = lines(i).point1;
        lines_expand(count).point2 = lines(i).point2;
        lines_expand(count).lineclass = 3;
    end
    if lines(i).lineclass1==0 && lines(i).lineclass2==0 && lines(i).lineclass3==0
        count = count + 1;
        lines_expand(count).point1 = lines(i).point1;
        lines_expand(count).point2 = lines(i).point2;
        lines_expand(count).lineclass = 0;
    end
end

lines_expand(count+1:end) = [];