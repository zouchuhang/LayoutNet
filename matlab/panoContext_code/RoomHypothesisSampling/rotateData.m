function [ oData ] = rotateData( iData )
%ROTATEDATA Rotate a structure
%   Detailed explanation goes here
numData = length(iData);
oData = repmat(iData(1), [numData 1]);
for did = 1:numData
    %% vp
    if isfield(iData, 'vp') && ~isempty(iData.vp)
        vps = iData(did).vp;
        selvp = vps(3:-1:1,:);
        R = diag([1 1 1])/(selvp'); %R * old coords = new coords;
        oData(did).vp = rotatePoint(vps, R);
        oData(did).vp(abs(oData(did).vp)<10^-5) = 0;
    end
    %% lines
%     lines = iData(did).lines;
    if isfield(iData, 'lines') && ~isempty(iData.lines)
        oData(did).lines = rotateLines( iData(did).lines, R );
        if isfield(iData(1), 'denseLines') && ~isempty(iData.denseLines)
            oData(did).denseLines = rotateLines( iData(did).denseLines, R);
        end
    end
    %% hypGT
%     hypGT = iData(did).hypGT;
    if isfield(iData, 'hypGT') && ~isempty(iData.hypGT)
        oData(did).hypGT = rotateHyps( iData(did).hypGT, R);
    end
    %% hyps
    if isfield(iData, 'hyps') && ~isempty(iData.hyps)
        oData(did).hyps = rotateHyps( iData(did).hyps, R);
    end
    
end
end

