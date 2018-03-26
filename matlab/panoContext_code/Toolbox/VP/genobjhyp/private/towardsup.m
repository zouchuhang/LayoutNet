function p = towardsup(p0, vp, amount)
dir = vp{1} - p0;
dir = dir / norm(dir);
if dir(2) > 0, dir = -dir; end
p = p0 + amount * dir;
