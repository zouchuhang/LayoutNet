function lr = determinepointlr(lc, p1, p2, vp)
% lr = 1 if p2 is on left of p1 (left means left, up, or front)
% lr = 2 if p2 is on right of p1 (right means right, down, or back)

if lc == 1
    if p1(2) < p2(2)
        lr = 2;
%         pl = p1;  pr = p2;
    else
        lr = 1;
%         pl = p2;  pr = p1;
    end
elseif lc == 2
    if p1(1) < p2(1)
        lr = 2;
%         pl = p1; pr = p2;
    else
        lr = 1;
%         pl = p2;  pr = p1;
    end
elseif lc == 3
    if acos( (p1-vp{3}) * (p2-vp{3})' ) > pi/2
        lr = 0;
    elseif norm(vp{3} - p1) > norm(vp{3} - p2)
        lr = 2;
%         pl = p1; pr = p2;
    else
        lr = 1;
%         pl = p2;  pr = p1;
    end
end

