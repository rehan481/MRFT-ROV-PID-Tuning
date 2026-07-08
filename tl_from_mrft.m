function [tl, Kp, Ki, Kd] = tl_from_mrft(beta_exp_vec, amp_exp_vec, freqHz_exp_vec, h, mode)
% TYREUS–LUYBEN (TL) GAINS FROM MRFT DATA
%   [tl, Kp, Ki, Kd] = tl_from_mrft(beta_exp_vec, amp_exp_vec, freqHz_exp_vec, h, mode)
%
% Uses MRFT measurements near beta≈0 (standard relay) to recover the
% ultimate gain/period and then applies Tyreus–Luyben tuning rules.
% Returns a struct 'tl' with P/PI/PID gains and convenience outputs
% [Kp,Ki,Kd] for the requested 'mode' ('P'|'PI'|'PID', default 'PID').
%
% Inputs:
%   beta_exp_vec     : vector of tested beta values
%   amp_exp_vec      : corresponding |e(t)| amplitudes at steady MRFT limit cycle
%   freqHz_exp_vec   : corresponding oscillation frequencies (Hz)
%   h                : MRFT relay amplitude
%   mode             : 'P' | 'PI' | 'PID'  (optional, default 'PID')
%
% Notes:
%   - Ultimate gain/period estimated as:
%       Ku = 4*h/(pi*A),  Tu = 1/fu,
%     using the MRFT point closest to beta=0.
%   - Tyreus–Luyben rules (parallel/ideal PID form):
%       PI : Kp = Ku/3.2,  Ti = 2.2*Tu       -> Ki = Kp/Ti, Kd = 0
%       PID: Kp = Ku/2.2,  Ti = 2.2*Tu, Td = Tu/6.3
%            -> Ki = Kp/Ti, Kd = Kp*Td
%       P  : Kp = Ku/2.2,  Ki = 0,           Kd = 0   (P-only TL variant)
%
% Example:
%   [tl, Kp,Ki,Kd] = tl_from_mrft(beta_vec, A_vec, f_vec, mrft_h, 'PID');
%   fprintf('TL PID: Kp=%.3g, Ki=%.3g, Kd=%.3g\n', Kp,Ki,Kd);

    if nargin < 5 || isempty(mode), mode = 'PID'; end
    assert(isequal(numel(beta_exp_vec), numel(amp_exp_vec), numel(freqHz_exp_vec)), ...
        'beta/amp/freq vectors must have equal length');

    % --- Find the MRFT point nearest to beta = 0 (standard relay case) ---
    [~, idx0] = min(abs(beta_exp_vec));   % exact 0 if present, else nearest
    beta0 = beta_exp_vec(idx0);
    A     = amp_exp_vec(idx0);            % error amplitude at oscillation
    fu    = freqHz_exp_vec(idx0);         % Hz
    Tu    = 1/max(eps, fu);               % ultimate period (s), guard fu=0
    Ku    = 4*h/(pi*max(eps, A));         % ultimate gain from DF (guard A=0)

    % --- Tyreus–Luyben rules ---
    % Precompute PI and PID constants
    % PI
    Kp_PI = Ku/3.2;
    Ti_PI = 2.2*Tu;
    Ki_PI = Kp_PI / Ti_PI;
    % PID
    Kp_PID = Ku/2.2;
    Ti_PID = 2.2*Tu;
    Td_PID = Tu/6.3;
    Ki_PID = Kp_PID / Ti_PID;
    Kd_PID = Kp_PID * Td_PID;
    % P (conservative, aligned with TL PID Kp scaling)
    Kp_P = Ku/2.2;

    % Package all options
    tl = struct();
    tl.beta_used = beta0;
    tl.A = A; tl.fu = fu; tl.Tu = Tu; tl.Ku = Ku;

    tl.P   = struct('Kp', Kp_P,   'Ki', 0,       'Kd', 0);
    tl.PI  = struct('Kp', Kp_PI,  'Ki', Ki_PI,   'Kd', 0);
    tl.PID = struct('Kp', Kp_PID, 'Ki', Ki_PID,  'Kd', Kd_PID);

    % Convenience return for selected mode
    switch upper(string(mode))
        case "P"
            Kp = tl.P.Kp;   Ki = 0;         Kd = 0;
        case "PI"
            Kp = tl.PI.Kp;  Ki = tl.PI.Ki;  Kd = 0;
        otherwise % 'PID'
            Kp = tl.PID.Kp; Ki = tl.PID.Ki; Kd = tl.PID.Kd;
    end
end
