function [phat, pci] = mlecustom(data,varargin)
%MLE Maximum likelihood estimation for custom univariate distributions.
%
%   See help for MLE.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.7 $  $Date: 2004/04/04 03:42:19 $

% Process any optional input arguments.
pnames = {'pdf' 'cdf' 'logpdf' 'logsf' 'nloglf' 'start' 'censoring' 'frequency' ...
          'alpha' 'options' 'lowerbound' 'upperbound' 'optimfun' 'optimoptions'};
defaults =  {[] [] [] [] [] [] zeros(size(data)) ones(size(data)) ...
             0.05 [] [] [] 'fminsearch' []};
[eid,errmsg,pdf,cdf,logpdf,logsf,nloglf,start,cens,freq,alpha,userOpts,...
    lb,ub,optimFun,optimOpts] = statgetargs(pnames, defaults, varargin{:});
if ~isempty(eid)
   error(sprintf('stats:mle:%s',eid),errmsg);
end

userOpts = statset(statset('mlecustom'), userOpts);
haveGrad = strcmp(userOpts.GradObj, 'on');
checkFunVals = strcmp(userOpts.FunValCheck, 'on');
delta = userOpts.DerivStep;

% Determine how the distribution is specified, and get handles to the
% specified functions.
if ~isempty(nloglf)
    if isa(nloglf,'function_handle')
        nloglfFun = nloglf;
        nloglfAddArgs = {};
    elseif iscell(nloglf) && isa(nloglf{1},'function_handle')
        nloglfFun = nloglf{1};
        nloglfAddArgs = nloglf(2:end);
    else
        error('stats:mle:InvalidNloglf','The ''nloglf'' parameter value must be a function handle or a cell array\ncontaining a function handle.');
    end
    funType = 'nloglf';
    % Assume we'll use llf_nloglf.  Will switch to using llf_nloglf_diff if
    % optimFun is 'fmincon', haveGrad is false, and optimOpts.GradObj is 'on'.
    llf = @llf_nloglf;
    
    funArgs = {nloglfFun nloglfAddArgs{:}};

elseif ~isempty(logpdf)
    if isa(logpdf,'function_handle')
        logpdfFun = logpdf;
        logpdfAddArgs = {};
    elseif iscell(logpdf) && isa(logpdf{1},'function_handle')
        logpdfFun = logpdf{1};
        logpdfAddArgs = logpdf(2:end);
    else
        error('stats:mle:InvalidLogpdf','The ''logpdf'' parameter value must be a function handle or a cell array\ncontaining a function handle.');
    end
    if ~isempty(logsf)
        if isa(logsf,'function_handle')
            logsfFun = logsf;
            logsfAddArgs = {};
        elseif iscell(logsf) && isa(logsf{1},'function_handle')
            logsfFun = logsf{1};
            logsfAddArgs = logsf(2:end);
        else
            error('stats:mle:InvalidLogsf','The ''logsf'' parameter value must be a function handle or a cell array\ncontaining a function handle.');
        end
        funType = 'logpdflogsf';
    elseif isempty(logsf) && sum(cens) == 0
        logsfFun = [];
        logsfAddArgs = {};
        funType = 'logpdf';
    else
        error('stats:mle:LogSFRequired','You must provide a log SF along with the log PDF if the data include censoring.');
    end
    llf = @llf_logpdflogsf;
    haveGrad = false;
    
    fun1Args = {logpdfFun logpdfAddArgs{:}};
    fun2Args = {logsfFun logsfAddArgs{:}};

elseif ~isempty(pdf)
    if isa(pdf,'function_handle')
        pdfFun = pdf;
        pdfAddArgs = {};
    elseif iscell(pdf) && isa(pdf{1},'function_handle')
        pdfFun = pdf{1};
        pdfAddArgs = pdf(2:end);
    else
        error('stats:mle:InvalidPdf','The ''pdf'' parameter value must be a function handle or a cell array\ncontaining a function handle.');
    end
    if ~isempty(cdf)
        if isa(cdf,'function_handle')
            cdfFun = cdf;
            cdfAddArgs = {};
        elseif iscell(cdf) && isa(cdf{1},'function_handle')
            cdfFun = cdf{1};
            cdfAddArgs = cdf(2:end);
        else
            error('stats:mle:InvalidCdf','The ''cdf'' parameter value must be a function handle or a cell array\ncontaining a function handle.');
        end
        funType = 'pdfcdf';
    elseif isempty(cdf) && sum(cens) == 0
        cdfFun = [];
        cdfAddArgs = {};
        funType = 'pdf';
    else
        error('stats:mle:CDFRequired','You must provide a CDF with the PDF if the data include censoring.');
    end
    llf = @llf_pdfcdf;
    haveGrad = false;

    fun1Args = {pdfFun pdfAddArgs{:}};
    fun2Args = {cdfFun cdfAddArgs{:}};
    
else
    error('stats:mle:DistFunsRequired','You must provide function handles as values for either the ''pdf'' and ''cdf'' parameters, or the ''logpdf''\nand ''logsf'' parameters, or the ''nloglf'' parameter.');
end

% Determine the size of the parameter vector from the (required) specified
% initial parameter values.
if ~isempty(start)
    nparams = numel(start);
else
    error('stats:mle:StartRequired','You must supply the ''start'' parameter value.');
end

% Make sure specified parameter bounds have correct sizes.  The defaults
% are set here because they depend on the size of the params.
if isempty(lb)
    lb = repmat(-Inf, size(start));
elseif ~isequal(size(start), size(lb))
    error('stats:mle:BndsSizeMismatch','The values of the ''lowerbound'' and ''start'' parameters must be the same size.');
end
if isempty(ub)
    ub = repmat(Inf, size(start));
elseif ~isequal(size(start), size(ub))
    error('stats:mle:BndsSizeMismatch','The values of the ''upperbound'' and ''start'' parameters must be the same size.');
end
if ~all((lb <= start) & (start <= ub))
    error('stats:mle:StartOutOfRange','The value of the ''start'' parameter must satisfy the upper and lower parameter bounds.');
end

% Now that start is valid, check once for errors thrown by the
% user-supplied functions.
if isequal(funType, 'nloglf')
    checkFunErrs('nloglf',nloglfFun,start,data,cens,freq,nloglfAddArgs,haveGrad);
    
else
    % For the PDF/CDF or logPDF/logSF forms, divide the data up into censored
    % and uncensored observations to save work later.
    c = (cens == 0);
    uncensData = data(c);
    uncensFreq = freq(c);
    censData = data(~c);
    censFreq = freq(~c);
    
    switch funType
    case {'logpdflogsf' 'logpdf'}
        checkFunErrs('logpdf',logpdfFun,start,uncensData,[],[],logpdfAddArgs);
        checkFunErrs('logsf',logsfFun,start,censData,[],[],logsfAddArgs);
    otherwise
        checkFunErrs('pdf',pdfFun,start,uncensData,[],[],pdfAddArgs);
        checkFunErrs('cdf',cdfFun,start,censData,[],[],cdfAddArgs);
    end
end

% Use fminsearch (the default), with any specified options.
if isequal(optimFun, 'fminsearch')
    opts = optimset(userOpts, 'GradObj','off', 'FunValCheck','off');
    if ~isempty(optimOpts), opts = optimset(opts, optimOpts); end
    switch funType
    case 'nloglf'
        % Call llf_nloglf.  fminsearch will never ask for a gradient.
        [phat,nll,err,output] = ...
            fminsearch(llf,start,opts,data,cens,freq,funArgs,checkFunVals,lb,ub);
    otherwise
        % Call either llf_pdfcdf or llf_logpdflogsf.  Neither can return a
        % gradient, and fminsearch will never ask for one.
        [phat,nll,err,output] = ...
            fminsearch(llf,start,opts,uncensData,censData,uncensFreq,censFreq,fun1Args,fun2Args,checkFunVals,lb,ub);
    end

% If requested, use fmincon, with any specified options.
elseif isequal(optimFun, 'fmincon')
    opts = optimset(userOpts,'LargeScale','on', 'GradObj','on', 'PrecondBandWidth',Inf, 'FunValCheck','off');
    if ~isempty(optimOpts), opts = optimset(opts, optimOpts); end
    
    % Use TolBnd to approximate open bounds.
    lb = lb + userOpts.TolBnd;
    ub = ub - userOpts.TolBnd;

    switch funType
    case 'nloglf'
        % If nloglf can return a gradient, or if fmincon does not expect
        % one, call llf_nloglf.
        if haveGrad || strcmp(opts.GradObj, 'off')
            [phat,nll,err,output] = ...
                fmincon(llf,start,[],[],[],[],lb,ub,[],opts,data,cens,freq,funArgs,checkFunVals);
        
        % If nloglf cannot return a gradient, and fmincon expects one, call
        % llf_nloglf_diff to compute a finite difference approximation.
        % This is the default (large scale) situation for nloglf/fmincon.
        else
            llf = @llf_nloglf_diff;
            [phat,nll,err,output] = ...
                fmincon(llf,start,[],[],[],[],lb,ub,[],opts,data,cens,freq,funArgs,checkFunVals,delta);
        end
        
    otherwise
        if strcmp(opts.GradObj, 'on')
            % Call a wrapper that will in turn call either llf_pdfcdf or
            % llf_logpdflogsf.  Neither can return a gradient, so, the wrapper
            % will have do it numerically.  This is the default (large scale)
            % situation for pdf/cdf/fmincon and logpdf/logsf/fmincon.
            [phat,nll,err,output] = ...
                fmincon(@llf_diff,start,[],[],[],[],lb,ub,[],opts,uncensData,censData,uncensFreq,censFreq,llf,fun1Args,fun2Args,checkFunVals,delta);
        else
            % Call either llf_pdfcdf or llf_logpdflogsf directly.  Neither can
            % return a gradient, but because opts.GradObj is 'off', fmincon
            % will never ask for one.
            [phat,nll,err,output] = ...
                fmincon(llf,start,[],[],[],[],lb,ub,[],opts,uncensData,censData,uncensFreq,censFreq,fun1Args,fun2Args,checkFunVals);
        end
    end

else
    error('stats:mle:IllegalOptimFunValue','The value of the ''optimfun'' parameter must be ''fminsearch'' or ''fmincon''.');
end

if (err == 0)
    % the optimizer may print its own output text; in any case give something
    % more statistical here, controllable via warning IDs.
    if output.funcCount >= opts.MaxFunEvals
        wmsg = 'Maximum likelihood estimation did not converge.  Function evaluation limit exceeded.';
    else
        wmsg = 'Maximum likelihood estimation did not converge.  Iteration limit exceeded.';
    end
    warning('stats:mle:IterOrEvalLimit',wmsg);
elseif err < 0
    error('stats:mle:MLEsNotFound','Unable to reach a maximum likelihood solution.');
end

% Compute the hessian with respect to the parameters, and invert it to get
% an estimate of the asymptotic covariance matrix of the param estimates.
if nargout > 1
    % The default FD step size is different for mlecov than for mlecustom,
    % because the former is computing a hessian, not a gradient.  We won't
    % tell mlecov what to use.
    userOpts.DerivStep = [];
    switch funType
    case 'nloglf'
        acov = mlecov(phat, data, 'nloglf',nloglf, 'cens',cens, 'freq',freq, 'options',userOpts);
    case {'logpdflogsf' 'logpdf'}
        acov = mlecov(phat, data, 'logpdf',logpdf, 'logsf',logsf, 'cens',cens, 'freq',freq, 'options',userOpts);
    case {'pdfcdf' 'pdf'}
        acov = mlecov(phat, data, 'pdf',pdf, 'cdf',cdf, 'cens',cens, 'freq',freq, 'options',userOpts);
    end

    % Compute CIs using a normal approximation for phat.
    probs = [alpha/2; 1-alpha/2];
    se = sqrt(diag(acov))';
    pci = norminv(repmat(probs,1,numel(phat)), [phat; phat], [se; se]);
end


%==========================================================================

function [nll,ngrad] = llf_nloglf(params, data, cens, freq, nloglfArgs, checkVals, lb, ub)
% Given a function handle to a negative log LF, evaluate the negative
% log-likelihood for PARAMS given DATA.

nloglfFun = nloglfArgs{1};
nloglfAddArgs = nloglfArgs(2:end);

% Bounds checking needed when function called by fminsearch.  When called
% by fmincon, the last two args will not be there.
if nargin > 6
    if any(params<=lb | ub<=params)
        nll = Inf;
        return
    end
end

if nargout == 1
    nll = feval(nloglfFun, params, data, cens, freq, nloglfAddArgs{:});
else
    [nll,ngrad] = feval(nloglfFun, params, data, cens, freq, nloglfAddArgs{:});
end

% Make sure returned llf values are valid.
if checkVals
    if ~isfinite(nll)
        error('stats:mle:NonfiniteNloglfVal','The NLOGLF function returned a NaN or infinite log-likelihood value.');
    end
    if nargout == 2
        if any(~isfinite(ngrad))
            error('stats:mle:NonfiniteNloglfGrad','The NLOGLF function returned NaN or infinite gradient values.');
        end
    end
end


%==========================================================================

function [nll,ngrad] = llf_nloglf_diff(params, data, cens, freq, nloglfArgs, checkVals, delta)
% Given a function handle to a negative log LF, evaluate the negative
% log-likelihood for PARAMS given DATA, and approximate its gradient using
% central differences.

nloglfFun = nloglfArgs{1};
nloglfAddArgs = nloglfArgs(2:end);

% No need to bounds check, this function only ever called by fmincon.

% Evaluate the log-likelihood itself, using the specified nlogLF.
nll = feval(nloglfFun, params, data, cens, freq, nloglfAddArgs{:});

% Approximate the gradient with central differences.
if nargout > 1
    deltaparams = delta*max(abs(params), 1); % limit smallest absolute step
    nparams = length(params);

    e = zeros(1,nparams);
    ngrad = zeros(size(params));
    for j=1:nparams
        e(j) = deltaparams(j);
        ngrad(j) = feval(nloglfFun, params+e, data, cens, freq, nloglfAddArgs{:}) ...
                 - feval(nloglfFun, params-e, data, cens, freq, nloglfAddArgs{:});
        e(j) = 0;
    end

    % Normalize by increment to get derivative estimates.
    ngrad = ngrad ./ (2 * deltaparams);
end

% Make sure returned llf values are valid.
if checkVals
    if ~isfinite(nll)
        error('stats:mle:NonfiniteNloglfVal','The NLOGLF function returned a NaN or infinite log-likelihood value.');
    end
    if nargout == 2
        if any(~isfinite(ngrad))
            error('stats:mle:NonfiniteNloglfVal','The NLOGLF function returned a NaN or infinite log-likelihood value.');
        end
    end
end


%==========================================================================

function nll = llf_logpdflogsf(params, uncensData, censData, uncensFreq, censFreq, logpdfArgs, logsfArgs, checkVals, lb, ub)
% Given function handles to a log PDF and a log SF, evaluate the negative
% log-likelihood for PARAMS given DATA.

logpdfFun = logpdfArgs{1};
logpdfAddArgs = logpdfArgs(2:end);
logsfFun = logsfArgs{1};
logsfAddArgs = logsfArgs(2:end);

% Bounds checking needed when function called by fminsearch.  When called
% by fmincon, the last two args will not be there.
if nargin > 8
    if any(params<=lb | ub<=params)
        nll = Inf;
        return
    end
end

% Log-likelihood = logPDF(uncensored values) + logSF(censored values)
%
% First, evaluate the specified logPDF of the uncensored data.
paramsCell = num2cell(params);
logpdfVals = feval(logpdfFun, uncensData, paramsCell{:}, logpdfAddArgs{:});

% Make sure returned logpdf values are valid.
if checkVals
    if any(~isfinite(logpdfVals))
        error('stats:mle:NonfiniteLogpdfVal','The LOGPDF function returned NaN or infinite values.');
    end
end

% Compute negative log-likelihood from uncensored values, using
% frequencies.
nll = -sum(uncensFreq.*logpdfVals);

% If there is censoring, evaluate the specified logSF of the censored data.
if ~isempty(censData)
    logsfVals = feval(logsfFun, censData, paramsCell{:}, logsfAddArgs{:});

    % Make sure returned logsf values are valid.
    if checkVals
        if any(~isfinite(logsfVals))
            error('stats:mle:NonfiniteLogsfVal','The LOGSF function returned NaN or infinite values.');
        elseif any(logsfVals > 0)
            error('stats:mle:PositiveLogsfVal','The LOGSF function returned positive values.');
        end
    end

    % Update negative log-likelihood with censored values, using
    % frequencies.
    nll = nll - sum(censFreq.*logsfVals);
end


%==========================================================================

function nll = llf_pdfcdf(params, uncensData, censData, uncensFreq, censFreq, pdfArgs, cdfArgs, checkVals, lb, ub)
% Given function handles to a PDF and a CDF, evaluate the negative
% log-likelihood for PARAMS given DATA.

pdfFun = pdfArgs{1};
pdfAddArgs = pdfArgs(2:end);
cdfFun = cdfArgs{1};
cdfAddArgs = cdfArgs(2:end);

% Bounds checking needed when function called by fminsearch.  When called
% by fmincon, the last two args will not be there.
if nargin > 8
    if any(params<=lb | ub<=params)
        nll = Inf;
        return
    end
end

% Log-likelihood = log(PDF(uncensored values) + log(1-CDF(censored values))
%
% First, evaluate the specified PDF of the uncensored data.
paramsCell = num2cell(params);
pdfVals = feval(pdfFun, uncensData, paramsCell{:}, pdfAddArgs{:});

% Make sure returned pdf values are valid.
if checkVals
    if any(~isfinite(pdfVals))
        error('stats:mle:NonfinitePdfVal','The PDF function returned NaN or infinite values.');
    elseif any(pdfVals <= 0)
        error('stats:mle:NonpositivePdfVal','The PDF function returned negative or zero values.');
    end
end

% Compute negative log-likelihood from uncensored values, using
% frequencies.
nll = -sum(uncensFreq.*log(pdfVals));

% If there is censoring, evaluate the specified CDF of the censored data.
if ~isempty(censData)
    cdfVals = feval(cdfFun, censData, paramsCell{:}, cdfAddArgs{:});

    % Make sure returned cdf values are valid.
    if checkVals
        if any(~isfinite(cdfVals))
            error('stats:mle:NonfiniteCdfVal','The CDF function returned NaN or infinite values.');
        elseif any(cdfVals < 0)
            error('stats:mle:NegativeCdfVal','The CDF function returned negative values.');
        elseif any(cdfVals >= 1)
            error('stats:mle:GTOneCdfVal','The CDF function returned values greater than or equal to 1.');
        end
    end

    % Update negative log-likelihood with censored values, using
    % frequencies.
    nll = nll - sum(censFreq.*log(1-cdfVals));
end


%==========================================================================

function [nll,ngrad] = llf_diff(params, uncensData, censData, uncensFreq, censFreq, llf, fun1Args, fun2Args, checkVals, delta)
% A wrapper around llf_logpdflogsf or llf_pdfcdf to evaluate the (negative)
% log-likelihood for PARAMS given DATA, and approximate its gradient using
% central differences.  LLF is a function handle to either llf_logpdflogsf
% or llf_pdfcdf, that in turn call either a PDF and CDF, or a logPDF and
% logSF, respectively.

% No need to bounds check, this function only ever called by fmincon.

% Evaluate the log-likelihood itself, using either the PDF and CDF, or the
% logPDF and logSF.
nll = feval(llf, params, uncensData, censData, uncensFreq, censFreq, fun1Args, fun2Args, checkVals);

% Approximate the gradient with central differences.
if nargout >= 2
    deltaparams = delta*max(abs(params), 1); % limit smallest absolute step
    nparams = length(params);

    e = zeros(1,nparams);
    ngrad = zeros(size(params));
    for j=1:nparams
        e(j) = deltaparams(j);
        ngrad(j) = feval(llf, params+e, uncensData, censData, uncensFreq, censFreq, fun1Args, fun2Args, checkVals) ...
                 - feval(llf, params-e, uncensData, censData, uncensFreq, censFreq, fun1Args, fun2Args, checkVals);
        e(j) = 0;
    end

    % Normalize by increment to get derivative estimates.
    ngrad = ngrad ./ (2 * deltaparams);
end

% No need to check values, that's done in llf in the above feval calls.


%==========================================================================

function checkFunErrs(type,fun,params,data,cens,freq,addArgs,haveGrad)
%CHECKFUNERRS Check for errors in evaluation of user-supplied function

if isempty(fun), return; end

try
    switch type
    case 'nloglf'
        if haveGrad
            [nll,ngrad] = feval(fun, params, data, cens, freq, addArgs{:});
        else
            nll = feval(fun, params, data, cens, freq, addArgs{:});
        end
    otherwise
        paramsCell = num2cell(params);
        vals = feval(fun, data, paramsCell{:}, addArgs{:});
    end
catch
    switch type
    case 'nloglf', errID = 'stats:mle:NloglfError';
    case 'logpdf', errID = 'stats:mle:LogpdfError';
    case 'logsf',  errID = 'stats:mle:LogsfError';
    case 'pdf',    errID = 'stats:mle:PdfError';
    case 'cdf',    errID = 'stats:mle:CdfError';
    end
    error(errID, ['The following error occurred while trying to evaluate\nthe ', ...
                  'user-supplied %s function ''%s'':\n\n%s'], type,func2str(fun),lasterr);
end
