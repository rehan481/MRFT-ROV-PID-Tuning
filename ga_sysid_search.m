function [thopt, fval, details] = ga_sysid_search(cost, th0, varargin)
% GA_SYSID_SEARCH  GA-based search for generic parameter vector th.
% Usage:
%   [thopt, fval] = ga_sysid_search(cost, th0)
%   [thopt, fval] = ga_sysid_search(cost, th0, 'lb',-100*abs(th0), 'ub',100*abs(th0))
%
% Inputs
%   cost : @(th)-> scalar cost
%   th0  : initial guess vector
%
% Name-Value (optional, minimal)
%   'lb'             : lower bounds (default: -1e3)
%   'ub'             : upper bounds (default:  +1e3)
%   'PopulationSize' : GA population (default: 60)
%   'MaxGenerations' : GA generations (default: 120)
%   'UseParallel'    : true/false (default: false)
%   'Polish'         : local fminsearch after GA (default: true)
%   'PolishOptions'  : optimset(...) for fminsearch (default basic)
%   'GAOptions'      : optimoptions('ga',...) to override defaults
%   'SeedSigma'      : std-dev for initial population scatter (default: 0.25*max(1,|th0|))

    if exist('ga','file') ~= 2
        error('ga_sysid_search: Global Optimization Toolbox (ga) is required.');
    end

    nvars = numel(th0);
    lb = -1e3*ones(1,nvars);
    ub =  1e3*ones(1,nvars);
    popSize = 60;
    maxGen  = 120;
    usePar  = false;
    doPolish = true;
    seedSigma = 0.25*max(1,abs(th0(:).')).';
    gaOptsUser = [];
    fmOpts = optimset('Display','off','MaxIter',400,'MaxFunEvals',5000);

    for k = 1:2:numel(varargin)
        name = lower(string(varargin{k}));
        val  = varargin{k+1};
        switch name
            case "lb",              lb = val(:).';
            case "ub",              ub = val(:).';
            case "populationsize",  popSize = val;
            case "maxgenerations",  maxGen  = val;
            case "useparallel",     usePar  = logical(val);
            case "polish",          doPolish = logical(val);
            case "polishoptions",   fmOpts = val;
            case "gaoptions",       gaOptsUser = val;
            case "seedsigma",       seedSigma = val(:).';
            otherwise, error('Unknown option "%s".', name);
        end
    end
    lb = lb(:).'; ub = ub(:).';
    if any(ub<=lb), error('Each element of ub must be > lb.'); end

    gaOpts = optimoptions('ga', ...
        'Display','iter', ...
        'PopulationSize', popSize, ...
        'MaxGenerations', maxGen, ...
        'UseParallel', usePar, ...
        'InitialPopulationMatrix', seed_population(th0, seedSigma, popSize, lb, ub), ...
        'FunctionTolerance', 1e-8, ...
        'ConstraintTolerance', 1e-6);

    if ~isempty(gaOptsUser)
        fns = fieldnames(gaOptsUser);
        for i = 1:numel(fns)
            try, gaOpts.(fns{i}) = gaOptsUser.(fns{i}); end %#ok<TRYNC>
        end
    end

    safeCost = @(x) local_safe_cost(cost, x);
    [x_ga, fval_ga, gaexit, gaoutput, population, scores] = ...
        ga(safeCost, nvars, [], [], [], [], lb, ub, [], gaOpts); %#ok<ASGLU>

    th_star = x_ga(:).'; 
    f_star  = fval_ga;

    if doPolish
        [x_local, f_local] = fminsearch(@(x) local_safe_cost(cost,x), th_star, fmOpts);
        if isfinite(f_local) && f_local < f_star
            th_star = x_local(:).';
            f_star  = f_local;
        end
    end

    thopt = th_star;
    fval  = f_star;
    details = struct('x_ga', x_ga, 'fval_ga', fval_ga, ...
                     'gaoutput', gaoutput, 'population', population, ...
                     'scores', scores);
end

% ---- helpers ----
function J = local_safe_cost(cost, x)
    th = x(:).';
    J = cost(th);
    if ~isfinite(J) || isnan(J), J = 1e300; end
end

function M = seed_population(th0, sigma, popSize, lb, ub)
    n = numel(th0);
    th0 = th0(:).';
    sigma = sigma(:).';
    if numel(sigma) ~= n, sigma = 0.25*max(1,abs(th0)); end
    M = zeros(popSize, n);
    M(1,:) = th0;
    for i = 2:popSize
        cand = th0 + sigma.*randn(1,n);
        cand = max(min(cand, ub), lb);
        M(i,:) = cand;
    end
end
