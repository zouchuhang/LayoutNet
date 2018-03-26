function vp = find_single_vp(lines, THRES_THETA)

NUMITER = 200;

% setup line samples and vote array
numlines = length(lines);
maxcomb = numlines*(numlines-1)/2;

if NUMITER > maxcomb
    % RANSAC
    votearray = zeros(NUMITER, 3); % [l m vote]
    for i = 1:NUMITER
        votearray(i, 1:2) = randsample(numlines,2);
    end
else
    % exhaustive
    votearray = zeros(maxcomb, 3); % [l m vote]
    c = 0;
    for i = 1:numlines
        for j = i+1:numlines
            c = c+1;
            votearray(c,1) = i;
            votearray(c,2) = j;
        end
    end
end

% test each vote array
for i = 1:size(votearray,1)
	% randomly select 2 lines
    ll = votearray(i, 1:2);

    vp = line_intersect_forvp(lines(ll(1)), lines(ll(2)));

	votearray(i,3) = getvote(vp, lines, THRES_THETA);
end

% the one with most vote
% n = find(votearray(:,3)==max(votearray(:,3)), 1, 'first');
[c maxidx] = max(votearray(:,3));

vp = line_intersect_forvp(lines(votearray(maxidx,1)),lines(votearray(maxidx,2)));

