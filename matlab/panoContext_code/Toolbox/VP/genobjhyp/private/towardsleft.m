function p = towardsleft(p0, vp, amount)
dir = vp{2} - p0;
dir = dir / norm(dir);
if dir(1) > 0, dir = -dir; end
p = p0 + amount * dir;
