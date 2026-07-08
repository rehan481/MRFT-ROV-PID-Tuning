%% === ZN gains from MRFT at beta ~ 0 (Åström–Hägglund) ===
% Uses: beta_exp_vec, amp_exp_vec (= max|e|), freqHz_exp_vec, mrft_h
% Returns a struct 'zn' and optional direct [Kp,Ki,Kd] for PID.

function [zn, Kp, Ki, Kd] = zn_from_mrft(beta_exp_vec, amp_exp_vec, freqHz_exp_vec, h, mode)
    if nargin < 5 || isempty(mode), mode = 'PID'; end   % 'P' | 'PI' | 'PID'
    assert(isequal(numel(beta_exp_vec),numel(amp_exp_vec),numel(freqHz_exp_vec)), ...
        'beta/amp/freq vectors must have equal length');

    % Find β = 0 (or the nearest available sample)
    [~, idx0] = min(abs(beta_exp_vec));   % exact 0 if present, else nearest
    beta0 = beta_exp_vec(idx0);
    A     = amp_exp_vec(idx0);            % amplitude of error e(t) at oscillation
    fu    = freqHz_exp_vec(idx0);         % Hz
    Tu    = 1/max(eps, fu);               % ultimate period
    Ku    = 4*h/(pi*A);                   % ultimate gain from describing function

    % Package results
    zn = struct();
    zn.beta_used = beta0;
    zn.A = A; zn.fu = fu; zn.Tu = Tu; zn.Ku = Ku;

    % ZN rules
    zn.P.Kp  = 0.5*Ku;

    zn.PI.Kp = 0.45*Ku;
    Ti_PI    = 0.83*Tu;
    zn.PI.Ki = zn.PI.Kp / Ti_PI; zn.PI.Kd = 0;

    zn.PID.Kp = 0.6*Ku;
    Ti_PID    = 0.5*Tu;
    Td_PID    = 0.125*Tu;
    zn.PID.Ki = zn.PID.Kp / Ti_PID;
    zn.PID.Kd = zn.PID.Kp * Td_PID;

    % Convenience return for chosen mode
    switch upper(mode)
        case 'P',   Kp = zn.P.Kp;   Ki = 0;           Kd = 0;
        case 'PI',  Kp = zn.PI.Kp;  Ki = zn.PI.Ki;    Kd = 0;
        otherwise,  Kp = zn.PID.Kp; Ki = zn.PID.Ki;   Kd = zn.PID.Kd; % 'PID'
    end
end
% 
