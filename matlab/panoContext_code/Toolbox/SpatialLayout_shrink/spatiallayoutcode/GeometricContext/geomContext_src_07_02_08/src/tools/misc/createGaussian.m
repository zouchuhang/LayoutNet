function G = createGaussian(x0,y0,sigmaX,sigmaY,gSize)
    radius=double(uint8((gSize-1)/2));
    G=zeros(gSize);
    for (x=-radius:radius)
        for (y=-radius:radius)
            G(x+radius+1,y+radius+1)=1/(2*pi*sigmaX*sigmaY)*exp(-1*((x-x0)^2)/(2*sigmaX^2) - (y-y0)^2/(2*sigmaY^2));
        end
    end
    G=G./sum(sum(G));