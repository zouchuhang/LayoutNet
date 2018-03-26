function [ inside, max_intrude_CT, max_intrude_ID ] = insideCone( ccwCones, vc, tol )
%INSIDECONE Check if vectors are in a cone, in 3D space
%   Cone is formed by ccwCones, vectors should be counter-clockwise viewing
%   from original point; vc is the vector for checking; tol is a tolerance
%   value, set tol=0 for exact judgement.
if ~exist('tol','var')
    tol = 0;
end

%% get outward normal of cone surface
numEdge = size(ccwCones, 1);
normal = zeros(numEdge, 3);
for i = 1:numEdge-1
    normal(i,:) = cross( ccwCones(i,:), ccwCones(i+1,:), 2);
end
normal(numEdge,:) = cross( ccwCones(numEdge,:), ccwCones(1,:), 2);

normal = normal./repmat( sqrt(sum(normal.^2, 2)), 1, 3);

%% negative to all outward normal
numVC = size(vc, 1);
inside = true(numVC, 1);
dotprods= zeros(numVC, numEdge);
for i = 1:numEdge
    dotprods(:,i) = dot( vc, repmat(normal(i,:), numVC, 1), 2);
    valid = dotprods(:,i) < cos((90-tol)*pi/180);
    inside = inside & valid;
end

[B,I] = max(dotprods,[],2);
max_intrude_CT = pi/2-acos(B);
max_intrude_ID = I;
end

