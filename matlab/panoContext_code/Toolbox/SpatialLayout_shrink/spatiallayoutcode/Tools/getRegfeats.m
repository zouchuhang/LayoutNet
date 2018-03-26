function [featsum]=getRegfeats(th1min,th1max,th2min,th2max,intImg)

%gives feat sum in the trapezoid 
featsum =intImg(th1max,th2max)+intImg(th1min,th2min)-...
intImg(th1min,th2max)-intImg(th1max,th2min);
    

return;