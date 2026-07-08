%% ================== Local helpers: case_a only ======================
function [Wmod_handle, th0] = create_caseA_sysid_model(w_init)
% case_a model only: integrator * 1st-order * 2nd-order + delay
% th = [k1, T_hydr, T_th, zeta, tau]

if nargin < 1 || isempty(w_init)
    Tseed = 0.2;
else
    wmed  = max(1, median(w_init(:)));
    Tseed = 1 / wmed;      % ~ one-radian time constant
end

Wmod_handle = @(th, w) caseA_sysid(th, w);
th0         = [1.0,  Tseed,  Tseed,  0.7,  0.02];
end

function W = caseA_sysid(th, w)
% th = [k1, T_hydr, T_th, zeta, tau]
assert(numel(th)==5, ...
    'caseA: th must have 5 elements [k1, T_hydr, T_th, zeta, tau].');

k1   = abs(th(1));
Th   = abs(th(2));
Tth  = abs(th(3));
zeta = abs(th(4));
tau  = max(0, th(5));

w  = w(:);
jw = 1j*w;

den_s    = jw;                                 % s
den_hydr = Th * jw + 1;                        % (T_hydr s + 1)
den_th   = (Tth^2) * (jw.^2) + 2*zeta*Tth*jw + 1;
% (T_th^2 s^2 + 2ζ T_th s + 1)

W = (k1 ./ (den_s .* den_hydr .* den_th)) .* exp(-1j*w*tau);
end