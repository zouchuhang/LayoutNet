function lines = assign_lineclass(lines, vp)

THRES_THETA = 10;
if length(vp)>=1 && ~isempty(vp{1})
    lineclass1 = line_belongto_vp(lines, vp{1}, THRES_THETA);
else
    lineclass1 = zeros(1, length(lines));
end
if length(vp)>=2 && ~isempty(vp{2})
    lineclass2 = line_belongto_vp(lines, vp{2}, THRES_THETA);
else
    lineclass2 = zeros(1, length(lines));
end
if length(vp)>=3 && ~isempty(vp{3})
    lineclass3 = line_belongto_vp(lines, vp{3}, THRES_THETA);
else
    lineclass3 = zeros(1, length(lines));
end


for i = 1:length(lines)
	lines(i).lineclass1 = lineclass1(i);
	lines(i).lineclass2 = lineclass2(i);
	lines(i).lineclass3 = lineclass3(i);
	
	if lineclass1(i) + lineclass2(i) + lineclass3(i) == 1
		lines(i).lineclass = 1*lineclass1(i) + 2*lineclass2(i) + 3*lineclass3(i);
	else
		lines(i).lineclass = 0;
	end
end

