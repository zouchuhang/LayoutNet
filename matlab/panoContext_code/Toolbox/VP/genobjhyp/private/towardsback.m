function p = towardsback(p0, vp, amount)
dir = vp{3} - p0;
dir = dir / norm(dir);
p = p0 + amount * dir;
