function [gopt, fval, details] = ga_pid_search(cost, g0, varargin)
% GA_PID_SEARCH  Genetic-algorithm search for PID gains.
% Usage:
%   [gopt, fval] = ga_pid_search(cost, g0)
%   [gopt, fval] = ga_pid_search(cost, g0, 'lb',[0 0 0], 'ub',[200 50 20])
%   [gopt, fval] = ga_pid_search(cost, g0, 'GAOptions', optimoptions('ga',...))
%   [gopt, fval] = ga_pid_search(cost, g0, 'Polish', true, 'PolishOptions', optimset(...))
%
% Inputs
%   cost  : @(g) -> scalar cost (finite). Your cost already uses abs(g).
%   g0    : initial guess [Kp Ki Kd] (used to seed the GA population)
%
% Name-Value options (all optional)
%   'lb'             : lower bounds (default: [0 0 0])
%   'ub'             : upper bounds (default: [1e3 1e3 1e3])
%   'GAOptions'      : options created by optimoptions('ga', ...)   (default set below)
%   'PopulationSize' : convenience override for GA pop size         (default: 60)
%   'MaxGenerations' : convenience override for GA generations      (default: 120)
%   'UseParallel'    : true/false for GA                            (default: false)
%   'Polish'         : do a quick local fminsearch after GA         (default: true)
%   'PolishOptions'  : options for fminsearch                       (default set below)
%   'SeedSigma'      : 1x3 std-dev used to scatter around g0        (default: 0.25*max(1,|g0|))
%
% Outputs
%   gopt    : best gains [Kp Ki Kd] (returned as nonnegative)
%   fval    : best cost value
%   details : struct with GA outputs (x_ga, fval_ga, gaoutput, population, scores, polished flag)

    if exist('ga','file') ~= 2
        error('ga_pid_search: Global Optimization Toolbox (ga) is required.');
    end

    % --- defaults ---
    nvars = numel(g0);
    assert(nvars==3, 'Expected g0 as [Kp Ki Kd].');
    lb = zeros(1,nvars);
    ub = 1e3*ones(1,nvars);
    popSize = 60;
    maxGen  = 120;
    usePar  = false;
    doPolish = true;
    seedSigma = 0.25*max(1,abs(g0));  % scatter around g0

    gaOptsUser = [];
    fmOpts = optimset('Display','off','MaxIter',400,'MaxFunEvals',5000);

    % --- parse NV pairs ---
    for k = 1:2:numel(varargin)
        name = lower(string(varargin{k}));
        val  = varargin{k+1};
        switch name
            case "lb",              lb = val(:).';
            case "ub",              ub = val(:).';
            case "gaoptions",       gaOptsUser = val;
            case "populationsize",  popSize = val;
            case "maxgenerations",  maxGen  = val;
            case "useparallel",     usePar  = logical(val);
            case "polish",          doPolish = logical(val);
            case "polishoptions",   fmOpts = val;
            case "seedsigma",       seedSigma = val(:).';
            otherwise, error('Unknown option "%s".', name);
        end
    end
    lb = lb(:).'; ub = ub(:).';
    assert(numel(lb)==nvars && numel(ub)==nvars, 'lb/ub must match length of g0.');
    if any(ub<=lb), error('Each element of ub must be > lb.'); end

    % --- GA options (defaults, can be overridden by GAOptions) ---
    gaOpts = optimoptions('ga', ...
        'Display','iter', ...
        'PopulationSize', popSize, ...
        'MaxGenerations', maxGen, ...
        'UseParallel', usePar, ...
        'InitialPopulationMatrix', seed_population(g0, seedSigma, popSize, lb, ub), ...
        'FunctionTolerance', 1e-8, ...
        'ConstraintTolerance', 1e-6);

    if ~isempty(gaOptsUser)
        % merge: user options take precedence
        fns = fieldnames(gaOptsUser);
        for i = 1:numel(fns)
            try
                gaOpts.(fns{i}) = gaOptsUser.(fns{i});
            catch
                % ignore unknown fields
            end
        end
    end

    % --- safe cost (avoid NaN/Inf killing GA) ---
    safeCost = @(x) local_safe_cost(cost, x);

    % --- run GA ---
    [x_ga, fval_ga, gaexit, gaoutput, population, scores] = ...
        ga(safeCost, nvars, [], [], [], [], lb, ub, [], gaOpts); %#ok<ASGLU>

    g_star = x_ga(:).'; 
    f_star = fval_ga;

    % --- optional local polish with fminsearch ---
    polished = false;
    if doPolish
        [x_local, f_local] = fminsearch(@(x) local_safe_cost(cost,x), g_star, fmOpts);
        if isfinite(f_local) && (f_local < f_star)
            g_star = x_local(:).';
            f_star = f_local;
            polished = true;
        end
    end

    % return nonnegative gains (your CLTools also uses abs(g), but this is nicer)
    gopt = abs(g_star);
    fval = f_star;

    % details
    details = struct('x_ga', x_ga, 'fval_ga', fval_ga, ...
                     'gaoutput', gaoutput, 'population', population, ...
                     'scores', scores, 'polished', polished);
end

% ---------- helpers ----------
function J = local_safe_cost(cost, x)
    % Ensure row vector; protect optimizer from NaN/Inf by large finite penalty
    g = x(:).';
    J = cost(g);
    if ~isfinite(J) || isnan(J)
        J = 1e300;
    end
end

function M = seed_population(g0, sigma, popSize, lb, ub)
    % Create an initial population centered at g0 with Gaussian spread (clipped to [lb,ub])
    n = numel(g0);
    M = zeros(popSize, n);
    M(1,:) = g0(:).';  % first row is exactly g0
    for i = 2:popSize
        cand = g0(:).' + sigma(:).'.*randn(1,n);
        cand = max(min(cand, ub), lb);
        M(i,:) = cand;
    end
end
