function [ omap, panoOmap ] = computePanoOmap( scene, edge, xyz )
%COMPUTEPANOOMAP compute orientation map
%   edge: line segments in separate views
%   xyz: vanishing point
%   OUTPUT:
%   omap: orientation map of separate views
%   panoOmap: project to panorama theta phi coordinate

uv = xyz2uvN(xyz,1);
omap = edge;
numViews = length(edge);
[H,W,~] = size(edge(1).img);
for i = 1:numViews
    edgeLst = edge(i).edgeLst;
    if size(edgeLst,1)==0
        omap(i).img = zeros(H,W,3);
        continue;
    end
    lines = struct('point1',[],'point2',[],'length',[]);
    lines(1,size(edgeLst,1)) = struct('point1',[],'point2',[],'length',[]);
    for j = 1:size(edgeLst,1)
        lines(j).point1 = edgeLst(j,1:2);
        lines(j).point2 = edgeLst(j,3:4);
        lines(j).length = sqrt(sum((edgeLst(j,1:2)-edgeLst(j,3:4)).^2));
    end
    
    x = edge(i).vx; y = edge(i).vy;
    fov = edge(i).fov;
    ANGx = uv(:,1); ANGy = uv(:,2);
    % compute the radius of ball
    [imH, imW, ~] = size(edge(i).img);
    R = (imW/2) / tan(fov/2);

    % im is the tangent plane, contacting with ball at [x0 y0 z0]
    x0 = R * cos(y) * sin(x);
    y0 = R * cos(y) * cos(x);
    z0 = R * sin(y);

    % plane function: x0(x-x0)+y0(y-y0)+z0(z-z0)=0
    % view line: x/alpha=y/belta=z/gamma
    % alpha=cos(phi)sin(theta);  belta=cos(phi)cos(theta);  gamma=sin(phi)
    alpha = cos(ANGy).*sin(ANGx);
    belta = cos(ANGy).*cos(ANGx);
    gamma = sin(ANGy);

    % solve for intersection of plane and viewing line: [x1 y1 z1]
    division = x0*alpha + y0*belta + z0*gamma;
    x1 = R*R*alpha./division;
    y1 = R*R*belta./division;
    z1 = R*R*gamma./division;

    % vector in plane: [x1-x0 y1-y0 z1-z0]
    % positive x vector: vecposX = [cos(x) -sin(x) 0]
    % positive y vector: vecposY = [x0 y0 z0] x vecposX
    vec = [x1-x0 y1-y0 z1-z0];
    vecposX = [cos(x) -sin(x) 0];
    deltaX = (vecposX*vec') / sqrt(vecposX*vecposX') + (imW+1)/2;
    vecposY = cross([x0 y0 z0], vecposX);
    deltaY = (vecposY*vec') / sqrt(vecposY*vecposY') + (imH+1)/2;
    
    vp{1,1} = [deltaX(1) deltaY(1)];
    vp{2,1} = [deltaX(2) deltaY(2)];
    vp{3,1} = [deltaX(3) deltaY(3)];
    
    lines_orig = lines; 
    [lines lines_ex] = taglinesvp(vp, lines_orig);
    [omapmore, OMAP_FACTOR] = compute_omap(lines, vp, [H W 3]);
    omap(i).img = double(omapmore);
    omap(i).lines_orig = lines_orig;
    omap(i).lines = lines;
    omap(i).vp = vp;
    linesImg = zeros(H, W, 3);
    for j = 1:length(lines)
        lineclass = lines(j).lineclass;
        if lineclass==0
            continue;
        end
        x = linspace( lines(j).point1(1)+1, lines(j).point2(1)+1, 1000);
        y = linspace( lines(j).point1(2)+1, lines(j).point2(2)+1, 1000);
        xx = max( min( round(x), W), 1);
        yy = max( min( round(y), H), 1);
        index = sub2ind( [H W], yy, xx);        
        linesImg(H*W*(lineclass-1)+index) = 1;
    end
    omap(i).linesImg = linesImg;
    
%     roomhyp = sample_roomhyp(1000, lines_ex, vp, [H W 3]);
%     omap(i).roomhyp = roomhyp;
    
%     [ bestHyp ] = evaluateRoomHyp( omap(i) );
%     disp_room(roomhyp(randsample(length(roomhyp),10)), scene(i).img, 1); % display some
    
%     cuboidhyp_omap = generate_cuboid_from_omap(omapmore, vp, OMAP_FACTOR);
%     disp_cubes(cuboidhyp_omap, scene(i).img, 1); % display all
%% original
%     img = edge(i).img;
%     [lines linesmore] = compute_lines(img);
%     if length(lines)<=3
%         omap(i).img = zeros(H,W,3);
%         continue;
%     end
%     
%     [vp f] = compute_vp(lines, [H W 3]);
%     lines_orig = lines; 
%     [lines lines_ex] = taglinesvp(vp, lines_orig);
%     linesmore_orig = linesmore; 
%     [linesmore linesmore_ex] = taglinesvp(vp, linesmore_orig);
%     [omapmore, OMAP_FACTOR] = compute_omap(lines, vp, [H W 3]);
%     omap(i).img = double(omapmore);

end

panoOmap = combineViews( omap, 2048, 1024 );

end

