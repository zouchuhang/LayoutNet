function entropy  = get_entropy(p)

p = p / sum(p);
entropy = -sum(p.*log(p));