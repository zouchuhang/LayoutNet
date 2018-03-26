function ls = sample_line(lines)

sample_rate = 1; % sample every 5 pixel on line

% build intermediate datastructure ls
for i = 1:length(lines)
	n_sample = ceil( norm(lines(i).point1-lines(i).point2) / sample_rate );
	ls(i).n_sample = n_sample;
	
	ls(i).sample = [ ...
		linspace(lines(i).point1(1), lines(i).point2(1), n_sample)' ...
		linspace(lines(i).point1(2), lines(i).point2(2), n_sample)' ];

	ls(i).lineclass = repmat(lines(i).lineclass, n_sample, 1);
end
