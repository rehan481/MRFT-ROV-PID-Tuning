function [Kp, Ki, Kd, c1, c2, c3] = PIDgains_GM_Beta(mrft_h, beta_exp_vec, amp_exp_vec, freqHz_exp_vec, varargin)
% Single-shot MRFT→PID (GM scalar). Defaults: beta=-0.2, GM=2, c2=Inf (Ki=0).

    % Defaults
    beta = -0.2; GM = 2; c2 = Inf;

    % Minimal name-value parsing
    for k = 1:2:numel(varargin)
        name = lower(char(varargin{k}));
        val  = varargin{k+1};
        switch name
            case 'beta', beta = val;
            case 'gm',   GM   = val;
            case 'c2',   c2   = val;
        end
    end

    % Interpolate MRFT measurements at |beta|
    b  = abs(beta);
    a0 = interp1(beta_exp_vec,  amp_exp_vec,    b, 'linear', 'extrap');
    f0 = interp1(beta_exp_vec,  freqHz_exp_vec, b, 'linear', 'extrap');   % Hz
    % Derived quantities
    x  = -beta / sqrt(1 - beta^2);    % x = 2*pi*c3 from beta relation
    c3 = x / (2*pi);
    c1 = 1 / (GM * sqrt(1 + x^2));    % gain-margin constraint (c2->Inf)

    % Controller gains (homogeneous rules)
    Kc = c1 * (4*mrft_h) / (pi*a0);
    if isinf(c2)
        Ki = 0;
    else
        Ti = c2 / f0;                 % since Ti = c2 * 2π / Ω0 and Ω0 = 2π f0
        Ki = Kc / Ti;                 % = Kc * f0 / c2
    end
    Td = c3 / f0;                     % Td = c3 * 2π / Ω0 = c3 / f0

    Kp = Kc;
    Kd = Kc * Td;
end
