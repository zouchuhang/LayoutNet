function type = getjunctype(cornerpt, line1, line2)
% junctype: 1--2
%           |  |
%           3--4
% HORZONRIGHT = 1;
% HORZONLEFT = -1;
% HORZONTOP = 1;
% HORZONBOT = -1;
% OTHER = 0;

if     line1.lineclass==1 && line2.lineclass==2
    hl = line2; vl = line1;
elseif line1.lineclass==2 && line2.lineclass==1
    hl = line1; vl = line2;
elseif line1.lineclass==1 && line2.lineclass==3
    hl = line2; vl = line1;
elseif line1.lineclass==3 && line2.lineclass==1
    hl = line1; vl = line2;
elseif line1.lineclass==2 && line2.lineclass==3
    hl = line1; vl = line2;
elseif line1.lineclass==3 && line2.lineclass==2
    hl = line2; vl = line1;
else
    error('jadsoifadpjdofadfa');
end

%%
TOL = 5;
if     cornerpt(1)-TOL <  hl.point1(1) && cornerpt(1)-TOL  < hl.point2(1)
    lr = 'left';
elseif cornerpt(1)+TOL >= hl.point1(1) && cornerpt(1)+TOL >= hl.point2(1)
    lr = 'right';
else
    type = 0;
    return;
end

if     cornerpt(2)-TOL <  vl.point1(2) && cornerpt(2)-TOL <  vl.point2(2)
    ud = 'up';
elseif cornerpt(2)+TOL >= vl.point1(2) && cornerpt(2)+TOL >= vl.point2(2)
    ud = 'down';
else
    type = 0;
    return;
end

if     strcmp(lr,'left') && strcmp(ud,'up')
    type = 1;
elseif strcmp(lr,'right') && strcmp(ud,'up')
    type = 2;
elseif strcmp(lr,'left') && strcmp(ud,'down')
    type = 3;
elseif strcmp(lr,'right') && strcmp(ud,'down')
    type = 4;
else
    error('bug in getjunctype.m');
end    

% % % % % % %%
% % % % % % % lr = is_line_leftorright_of_point(lines, point, vp); % 1 if right, -1 if left, 0 if neither
% % % % % % if     is_line_leftorright_of_point(hl, vl.point1, vp)==1 && ...
% % % % % %        is_line_leftorright_of_point(hl, vl.point2, vp)==1
% % % % % %     lr = HORZONRIGHT;
% % % % % % elseif is_line_leftorright_of_point(hl, vl.point1, vp)==-1 && ...
% % % % % %        is_line_leftorright_of_point(hl, vl.point2, vp)==-1
% % % % % %     lr = HORZONLEFT;
% % % % % % else
% % % % % %     lr = OTHER;
% % % % % % end
% % % % % % 
% % % % % % if     is_line_upordown_of_point(hl, vl.point1)==1 && ...
% % % % % %        is_line_upordown_of_point(hl, vl.point2)==1 % 1 if up, -1 if down, 0 if neither
% % % % % %     ud = HORZONTOP;
% % % % % % elseif is_line_upordown_of_point(hl, vl.point1)==-1 && ...
% % % % % %        is_line_upordown_of_point(hl, vl.point2)==-1
% % % % % %     ud = HORZONBOT;
% % % % % % else
% % % % % %     ud = OTHER;
% % % % % % end
% % % % % % 
% % % % % % if     lr == HORZONRIGHT && ud == HORZONTOP
% % % % % %     type = 1;
% % % % % % elseif lr == HORZONRIGHT && ud == HORZONBOT
% % % % % %     type = 3;
% % % % % % elseif lr == HORZONLEFT && ud == HORZONTOP
% % % % % %     type = 2;
% % % % % % elseif lr == HORZONLEFT && ud == HORZONBOT
% % % % % %     type = 4;
% % % % % % else
% % % % % %     type = 0;
% % % % % % end


