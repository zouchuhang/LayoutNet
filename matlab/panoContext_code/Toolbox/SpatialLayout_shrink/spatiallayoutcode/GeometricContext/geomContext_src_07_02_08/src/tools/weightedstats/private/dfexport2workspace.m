function dfexport2workspace()
% DFEXPORT2WORKSPACE Helper file for the Distribution Fitting tool
%    DFEXPORT2WORKSPACE gets the saved evaluated results and passes them to
%    export2wsdlg

%   Copyright 2003-2004 The MathWorks, Inc. 
%   $Revision: 1.1.6.2 $

results = dfgetset('evaluateResults');
export2wsdlg({'Save evaluate results to a MATLAB variable named:'}, ...
             {'evaluateresults'}, {results});
