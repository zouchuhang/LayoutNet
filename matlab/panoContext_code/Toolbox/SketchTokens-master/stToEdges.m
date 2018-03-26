function E = stToEdges( S, nms, suppress )
% Convert sketch tokens to edges.
%
% USAGE
%  E = stToEdges( S, [nms], [suppress] )
%
% INPUTS
%  S          - [h x w x (nTokens+1)] sketch token probability maps
%  nms        - [1] if true apply non-maximum suppression to edges
%  suppress   - [1] if true suppress boundary confidence in order to
%                   compensate padding effect
%
% OUTPUTS
%  E          - [h x w] edge probability map
%
% EXAMPLE
%
% See also stDetect
%
% Sketch Token Toolbox     V0.95
% Copyright 2013 Joseph Lim [lim@csail.mit.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License [see bsd.txt]

    % get default arguments
    if(nargin<2 || isempty(nms)),
        nms=1;
    end
    if(nargin<3 || isempty(suppress)),
        suppress=1;
    end

    % extract edge probabilities
    E = 1-S(:,:,end);

    % apply nms
    if nms
        O=edgeOrient(E,4);
        E=edgeNms(E,O,1); 
    end
    
    % apply boundary suppression
    if suppress
        S=5;        
        for s=1:S,
            E([s end-s+1],:,:)=E([s end-s+1],:,:)*(s-1)/S;
        end
        for s=1:S,
            E(:,[s end-s+1],:)=E(:,[s end-s+1],:)*(s-1)/S;
        end
    end

end

function E = edgeNms( E, O, r )
    % suppress locations where edge is stronger in orthogonal direction
    E1=imPad(E,r+1,'replicate');
    Dx=cos(O);
    Dy=sin(O);
    
    [ht,wd]=size(E1);
    [cs,rs]=meshgrid(r+2:wd-r-1,r+2:ht-r-1);
    for i=-r:r,
        if i==0
            continue; 
        end
        cs0=i*Dx+cs;
        dcs=cs0-floor(cs0);
        cs0=floor(cs0);
        
        rs0=i*Dy+rs;
        drs=rs0-floor(rs0);
        rs0=floor(rs0);
        
        E2 = (1-dcs).*(1-drs) .* E1(rs0+0+(cs0-1)*ht);
        E2 = E2 + dcs.*(1-drs) .* E1(rs0+0+(cs0-0)*ht);
        E2 = E2 + (1-dcs).*drs .* E1(rs0+1+(cs0-1)*ht);
        E2 = E2 + dcs.*drs .* E1(rs0+1+(cs0-0)*ht);
        E(E*1.01<E2) = 0;
    end
end

function O = edgeOrient( E, r )
    % compute very approximate orientation map from edge map
    E2=convTri(E,r);
    f=[-1 2 -1];
    
    Dx=conv2(E2,f,'same');
    Dy=conv2(E2,f','same');
    F=conv2(E2,[1 0 -1; 0 0 0; -1 0 1],'same')>0;
    Dy(F)=-Dy(F);
    O=mod(atan2(Dy,Dx),pi);
end
