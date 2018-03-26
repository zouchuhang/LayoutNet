function [ ohyps ] = rotateHyps( hyps, R )
%ROTATEHYPS Rotate room layout hypothesis
%   
numHyps = length(hyps);
ohyps = hyps; %repmat(hyps(1), [numHyps 1]);
for hid = 1:numHyps
    ohyps(hid).extLine = rotateLines(hyps(hid).extLine, R);
    ohyps(hid).srcLine = rotateLines(hyps(hid).srcLine, R);
    ohyps(hid).hCorner = rotatePoint(hyps(hid).hCorner, R);
end



end

