function [lines lines_ex] = ...
    taglinesvp(vp, lines)

lines = assign_lineclass(lines, vp);

lines = compute_line_attributes(lines, vp);

lines_ex = expand_ambiguous_lineclass(lines);

lines_ex = compute_line_attributes(lines_ex, vp);

