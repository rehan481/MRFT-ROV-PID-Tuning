function info = print_sysid_model(model_type, th, varargin)
% PRINT_SYSID_MODEL  Pretty-print the identified transfer function, with normalization.
% Usage (pass the SAME trailing args you used with create_sysid_model):
%   info = print_sysid_model('caseB',     th_hat);
%   info = print_sysid_model('general',   th_hat, nz, np);
%   info = print_sysid_model('nth_delay', th_hat, np);
%   info = print_sysid_model('case_a',  th_hat);

% Normalization:
%   We divide BOTH numerator and denominator by 10^k, where k is chosen so
%   the largest coefficient magnitude is ~O(1). This keeps H(s) IDENTICAL,
%   just nicer to read. (Denominator may no longer be monic after scaling.)
%
% Returns:
%   info.num, info.den        % unscaled coefficients
%   info.num_s, info.den_s    % scaled (printed) coefficients
%   info.delay                % delay (if any)
%   info.scale10, info.k10    % scale factor used: divide by 10^k10
%   info.str_expanded         % printed string

mdl = lower(strtrim(model_type));
th  = th(:).';  % row

% local helpers
fmt = '%.6g';                      % nice, compact
sig = 6;                           % significant digits for rounding
rnd = @(x) round(x, sig, 'significant');
polystr = @(c,var) poly2str_fmt(c, var, fmt);

switch mdl
    case 'general'
        assert(numel(varargin) >= 2, 'general: provide nz and np.');
        nz = varargin{1};
        np = varargin{2};
        n_num = nz + 1;
        assert(numel(th) == (n_num + np), 'general: th must be (nz+1)+np.');

        num = rnd(th(1:n_num));
        den = rnd([1, th(n_num+1:end)]);    % monic before scaling

        % ---- normalize BOTH numerator & denominator by 10^k ----
        [num_s, den_s, k10, scale10] = normalize_b10(num, den);

        % pretty print (NO tf metadata)
        fprintf('\nIdentified GENERAL rational model (nz=%d, np=%d):\n', nz, np);
        fprintf('  H(s) = ( %s ) / ( %s )\n', polystr(num_s,'s'), polystr(den_s,'s'));
        %             if k10 ~= 0
        %                 fprintf('  [coefficients divided by 1e%d for readability; H(s) unchanged]\n', k10);
        %             end

        info.num = num; info.den = den; info.delay = 0;
        info.num_s = num_s; info.den_s = den_s;
        info.k10 = k10; info.scale10 = scale10;
        info.str_expanded = sprintf('( %s ) / ( %s )', polystr(num_s,'s'), polystr(den_s,'s'));

    case 'case_b'
        % th = [m, c1, ks, ka, wn, zeta]
        assert(numel(th)==6, 'case B: th must have 6 elements.');
        m    = abs(th(1));  c1 = abs(th(2));  ks = abs(th(3));
        ka   = abs(th(4));  wn = abs(th(5));  z  = abs(th(6));

        num = rnd(ka*wn^2);                % scalar
        den1 = rnd([m, c1, ks]);           % m s^2 + c1 s + ks
        den2 = rnd([1, 2*z*wn, wn^2]);     % s^2 + 2 z wn s + wn^2
        den  = rnd(conv(den1, den2));

        [num_s, den_s, k10, scale10] = normalize_b10(num, den);

        fprintf('\nIdentified CASE B model:\n');
        fprintf('  H(s) = ( %s ) / ( %s )\n', polystr(num_s,'s'), polystr(den_s,'s'));
        %             if k10 ~= 0
        %                 fprintf('  [coefficients divided by 1e%d for readability; H(s) unchanged]\n', k10);
        %             end

        info.num = num; info.den = den; info.delay = 0;
        info.num_s = num_s; info.den_s = den_s;
        info.k10 = k10; info.scale10 = scale10;
        info.str_expanded = sprintf('( %s ) / ( %s )', polystr(num_s,'s'), polystr(den_s,'s'));
    case 'case_a'
        % th = [k1, Th, Tth, zeta, tau]
        assert(numel(th)==5, 'case_a: th must be [k1, Th, Tth, zeta, tau].');

        k1   = abs(th(1));
        Th   = abs(th(2));
        Tth  = abs(th(3));
        zeta = abs(th(4));
        tau  = max(0, th(5));

        % H(s) = k1 / [ s * (Th s + 1) * (Tth^2 s^2 + 2 ζ Tth s + 1) ] * e^{-s tau}
        num  = rnd(k1);                               % scalar
        den1 = [1 0];                                 % s
        den2 = [Th 1];                                % (Th s + 1)
        den3 = [Tth^2, 2*zeta*Tth, 1];                % (Tth^2 s^2 + 2 ζ Tth s + 1)
        den  = rnd(conv(conv(den1, den2), den3));

        [num_s, den_s, k10, scale10] = normalize_b10(num, den);

        fprintf('\nIdentified CASE A model (hydraulic + thermal + delay):\n');
        fprintf('  H(s) = ( %s ) / ( %s ) * exp(-s*%s)\n', ...
            polystr(num_s,'s'), polystr(den_s,'s'), num2str(rnd(tau), fmt));

        info.num = num; info.den = den; info.delay = tau;
        info.num_s = num_s; info.den_s = den_s;
        info.k10 = k10; info.scale10 = scale10;
        info.str_expanded = sprintf('( %s ) / ( %s ) * exp(-s*%s)', ...
            polystr(num_s,'s'), polystr(den_s,'s'), num2str(rnd(tau), fmt));
    case {'case_c','casec'}
        % th = [k, wb, zb, T, z, tau]  (normalized body)
        k    = th(1);                 % keep sign
        wb   = abs(th(2));  zb = abs(th(3));
        Tact = abs(th(4));  z  = abs(th(5));
        tau  = max(0, th(6));

        % H(s) = k / [(s^2 + 2*zb*wb*s + wb^2) * (T^2 s^2 + 2*z*T s + 1)] * e^{-s tau}
        num  = rnd(k);                                % scalar
        den1 = [1, 2*zb*wb, wb^2];
        den2 = [Tact^2, 2*z*Tact, 1];
        den  = rnd(conv(den1, den2));

        % normalize both numerator & denominator by 10^k (cosmetic only)
        [num_s, den_s, k10, scale10] = normalize_b10(num, den);

        fprintf('\nIdentified CASE C model (normalized body):\n');
        fprintf('  H(s) = ( %s ) / ( %s ) * exp(-s*%s)\n', ...
            polystr(num_s,'s'), polystr(den_s,'s'), num2str(rnd(tau), fmt));

        % export details just like the other branches
        info.num = num; info.den = den; info.delay = tau;
        info.num_s = num_s; info.den_s = den_s;
        info.k10 = k10; info.scale10 = scale10;
        info.str_expanded = sprintf('( %s ) / ( %s ) * exp(-s*%s)', ...
            polystr(num_s,'s'), polystr(den_s,'s'), num2str(rnd(tau), fmt));
    case 'case_d'
        % th = [k, wb, zb, T, z, tau, w_z]
        assert(numel(th)==7, 'case_d: th must be [k, wb, zb, T, z, tau, w_z].');

        k    = th(1);                 % keep sign
        wb   = abs(th(2));  zb = abs(th(3));
        Tact = abs(th(4));  z  = abs(th(5));
        tau  = max(0, th(6));
        wz   = max(1e-6, abs(th(7)));

        % Body and actuator same as case_c
        den1 = [1, 2*zb*wb, wb^2];
        den2 = [Tact^2, 2*z*Tact, 1];
        den  = rnd(conv(den1, den2));

        % Zero: (1 + s/w_z) => num_zero = [1/w_z, 1]
        num_zero = [1/wz, 1];
        num      = rnd(k * num_zero);     % k*(1 + s/w_z)

        % normalize both numerator & denominator by 10^k (cosmetic only)
        [num_s, den_s, k10, scale10] = normalize_b10(num, den);

        fprintf('\nIdentified CASE D model (normalized body + real zero):\n');
        fprintf('  H(s) = ( %s ) / ( %s ) * exp(-s*%s)\n', ...
            polystr(num_s,'s'), polystr(den_s,'s'), num2str(rnd(tau), fmt));

        info.num = num; info.den = den; info.delay = tau;
        info.num_s = num_s; info.den_s = den_s;
        info.k10 = k10; info.scale10 = scale10;
        info.str_expanded = sprintf('( %s ) / ( %s ) * exp(-s*%s)', ...
            polystr(num_s,'s'), polystr(den_s,'s'), num2str(rnd(tau), fmt));


    case 'nth_delay'
        % th = [K, a_{np-1}, ..., a0, tau]
        assert(numel(varargin) >= 1, 'nth_delay: provide np.');
        np = varargin{1};
        assert(numel(th) == 1 + np + 1, 'nth_delay: th must be [K, a..., tau].');

        K   = th(1);
        a   = th(2:end-1);
        tau = max(0, th(end));

        num = rnd(K);                      % scalar
        den = rnd([1, a(:).']);            % monic before scaling

        [num_s, den_s, k10, scale10] = normalize_b10(num, den);

        fprintf('\nIdentified Nth-order all-pole with delay (np=%d):\n', np);
        fprintf('  H(s) = ( %s ) / ( %s ) * exp(-s*%s)\n', polystr(num_s,'s'), polystr(den_s,'s'), num2str(rnd(tau), fmt));
        %             if k10 ~= 0
        %                 fprintf('  [coefficients divided by 1e%d for readability; H(s) unchanged]\n', k10);
        %             end

        info.num = num; info.den = den; info.delay = tau;
        info.num_s = num_s; info.den_s = den_s;
        info.k10 = k10; info.scale10 = scale10;
        info.str_expanded = sprintf('( %s ) / ( %s ) * exp(-s*%s)', polystr(num_s,'s'), polystr(den_s,'s'), num2str(rnd(tau), fmt));

    otherwise

        error('Unknown model_type. Use ''caseB'', ''general'', ''nth_delay'', ''case_a'', or ''case_c''.');
end
fprintf('\n');
end

% ---------- helpers ----------
function [num_s, den_s, k10, scale10] = normalize_b10(num, den)
% Choose 10^k so max(|coef|)/10^k is in [1,10)
coeffs = [num(:); den(:)];
M = max(abs(coeffs));
if M==0
    k10 = 0; scale10 = 1;
else
    k10 = floor(log10(M));      % integer power
    scale10 = 10^k10;
end
num_s = num / scale10;
den_s = den / scale10;
end

function s = poly2str_fmt(c, var, fmt)
% Build polynomial string with custom numeric format
c = c(:).';
n = numel(c)-1;
parts = {};
for i = 1:numel(c)
    a = c(i);
    p = n-(i-1);
    if abs(a) < eps, continue; end
    % sign handling
    signStr = '';
    if isempty(parts)
        % first term: keep sign on the number itself
        coefStr = num2str(a, fmt);
    else
        if a >= 0
            signStr = ' + ';
            coefStr = num2str(abs(a), fmt);
        else
            signStr = ' - ';
            coefStr = num2str(abs(a), fmt);
        end
    end
    % variable power
    if p > 1
        term = sprintf('%s%s %s^%d', signStr, coefStr, var, p);
    elseif p == 1
        term = sprintf('%s%s %s', signStr, coefStr, var);
    else
        term = sprintf('%s%s', signStr, coefStr);
    end
    parts{end+1} = term; %#ok<AGROW>
end
if isempty(parts)
    s = '0';
else
    s = strjoin(parts, '');
end
end
