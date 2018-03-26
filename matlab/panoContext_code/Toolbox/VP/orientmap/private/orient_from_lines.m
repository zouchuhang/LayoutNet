function [lineextimg c] = orient_from_lines(lines, vp, imgwidth, imgheight)


%%
ls = sample_line(lines);
linesamples = cat(1,ls(:).sample);
linesampleclass = cat(1,ls(:).lineclass);

%%
% lineextimg = cell(3,3);
% for i=1:9, lineextimg{i} = zeros(imgheight,imgwidth); end
lineextimg = cell(3,3,2);
for i=1:18, lineextimg{i} = zeros(imgheight,imgwidth); end

%%
% poly = extend_line(line, vp{1}, stoppinglines_sample, imgwidth, imgheight);
for i = 1:length(lines)
    lc = lines(i).lineclass;
    if lc~=0
    for extdir = setdiff(1:3, lc)
        targetdir = setdiff(1:3, [lc extdir]);
    
%         poly = extend_line_old(lines(i), vp{extdir}, linesamples(linesampleclass==targetdir,:), imgwidth, imgheight);
%         lineextimg{lc,extdir} = lineextimg{lc,extdir} + poly2mask(poly(:,1), poly(:,2), imgheight, imgwidth);
        
        poly = extend_line(lines(i), vp{extdir}, 1, linesamples(linesampleclass==targetdir,:), imgwidth, imgheight);
        lineextimg{lc,extdir,1} = lineextimg{lc,extdir,1} + poly2mask(poly(:,1), poly(:,2), imgheight, imgwidth);
        
        poly = extend_line(lines(i), vp{extdir}, -1, linesamples(linesampleclass==targetdir,:), imgwidth, imgheight);
        lineextimg{lc,extdir,2} = lineextimg{lc,extdir,2} + poly2mask(poly(:,1), poly(:,2), imgheight, imgwidth);
    end
    end
end


%%
function poly = extend_line(line, vp, towards_or_away, stoppinglines_sample, imgwidth, imgheight)
% towards_or_away: 1 or -1

p1 = line.point1;
p2 = line.point2;

curp1 = p1; curp2 = p2;
move_amount = 128;
while move_amount>=1
    [newp1 newp2 atvp] = move_line_towards_vp(curp1, curp2, vp, towards_or_away * move_amount);
    
    failed = 0;
    if atvp==1
%         move_amount = 0; % exit now.
        failed = 1;
        
    elseif (newp1(1)>imgwidth || newp1(1)<1 || newp1(2)>imgheight || newp1(2)<1) && ...
       (newp2(1)>imgwidth || newp2(1)<1 || newp2(2)>imgheight || newp2(2)<1)
        failed = 1;
        
    else
        isstop = inpolygon(stoppinglines_sample(:,1), stoppinglines_sample(:,2), ...
            [p1(1) p2(1) newp2(1) newp1(1) p1(1)], [p1(2) p2(2) newp2(2) newp1(2) p1(2)]);
        
        if any(isstop)
            failed = 1;
        end
    end
    
    if failed
        move_amount = move_amount/2;
    else
        curp1 = newp1;
        curp2 = newp2;
    end
end
% poly = [curp1(:)'; curp2(:)'];
poly = [p1(:)'; p2(:)'; curp2(:)'; curp1(:)'];

% curp1 = p1; curp2 = p2;
% move_amount = 32;
% while move_amount>=1
%     [newp1 newp2 atvp] = move_line_towards_vp(curp1, curp2, vp, -move_amount);
%     
%     failed = 0;
%     if atvp==1
%         move_amount = 0; % exit now.
%         
%     elseif (newp1(1)>imgwidth || newp1(1)<1 || newp1(2)>imgheight || newp1(2)<1) && ...
%        (newp2(1)>imgwidth || newp2(1)<1 || newp2(2)>imgheight || newp2(2)<1)
%         failed = 1;
%         
%     else
%         isstop = inpolygon(stoppinglines_sample(:,1), stoppinglines_sample(:,2), ...
%             [p1(1) p2(1) newp2(1) newp1(1) p1(1)], [p1(2) p2(2) newp2(2) newp1(2) p1(2)]);
%         
%         if any(isstop)
%             failed = 1;
%         end
%     end
%     
%     if failed
%         move_amount = move_amount/2;
%     else
%         curp1 = newp1;
%         curp2 = newp2;
%     end
% end
% poly = [poly; curp2(:)'; curp1(:)'];
% poly = [poly; poly(1,:)];

%%

%%
function [newp1 newp2 atvp] = move_line_towards_vp(linep1, linep2, vp, amount)

% d = dist_line_to_point(linep1, linep2, vp);
% r = amount / d;
n1 = norm(vp-linep1);
n2 = norm(vp-linep2);
dir1 = (vp - linep1) / n1;
dir2 = (vp - linep2) / n2;
ratio21 = n2 / n1;

% if n1>amount && n2<amount
%     fprintf('check');
% end

if n1 < amount
    newp1 = linep1;
    newp2 = linep2;
    atvp = 1;
else
    newp1 = linep1 + dir1 * amount;
    newp2 = linep2 + dir2 * amount * ratio21;
    atvp = 0;
end

%%
function [newp1 newp2 atvp] = move_line_towards_vp_old(linep1, linep2, vp, amount)

d = dist_line_to_point(linep1, linep2, vp);
r = amount / d;
if r > 1
    newp1 = vp;
    newp2 = vp;
    atvp = 1;
else
    newp1 = linep1 + (vp-linep1)*r;
    newp2 = linep2 + (vp-linep2)*r;
    atvp = 0;
end
