function displaySegmentGraph(im, segimage, adjlist)

stats = regionprops(segimage, 'Centroid');
centroids = cat(1, stats.Centroid);

segimage = double(segimage) / double(max(segimage(:)));

[gx, gy] = gradient(segimage);

edges = find((gx ~= 0) | (gy ~= 0));

[imh, imw, imb] = size(im);

for b = 1:imb
    im(edges+(b-1)*imh*imw) = 1*(b==1); % make edges red
end

hold off, imagesc(im), axis image, hold on

% convert adjacency matrix to adjacency list if necessary
if size(adjlist, 1) == size(adjlist, 2)
    adjmat = adjlist;
    [s1, s2] = find(adjmat);
    adjlist = [s1 s2];
    ind = find(adjlist(:, 1) < adjlist(:, 2));
    adjlist = adjlist(ind, :);
end

if size(adjlist, 2) == 2
    for k = 1:size(adjlist, 1)
        plot(centroids(adjlist(k, :), 1), centroids(adjlist(k, :), 2), 'g', 'LineWidth', 1);
    end
end
    


