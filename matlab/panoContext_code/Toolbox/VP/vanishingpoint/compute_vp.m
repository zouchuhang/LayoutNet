function [vp f] = compute_vp(lines, imgsize, camcalib)
% compute vanishing points 'vp' using 'lines'
% assign line class to 'lines' and 'linesmore'

vpmultitype = 0;

% if camcalib.provided==1     % true camera parameters provided
%     [normal_vec, vp] = vanishline(lines, camcalib);
% else                        % images from the web

    if vpmultitype==0
        [vp f] = vanish_from_minevidence(lines, imgsize);

    elseif vpmultitype==1
        [normal_vec, vp] = vanishline_weakorthoconst(lines, camcalib);

    elseif vpmultitype==2
        vpset = vanishline_multi(lines, camcalib);
        % script_dispmultivp(vpset, lines, linesmore, img);
        % figure;

    elseif vpmultitype==3
        vpset = vanishline_multi3(lines, camcalib);
        % script_dispmultivp(vpset, lines, linesmore, img);
        % figure;

    elseif vpmultitype==4
        vpset = vanishline_multi4(lines, camcalib);
        % script_dispmultivp(vpset, lines, linesmore, img);
    end

% end


