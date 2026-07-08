%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Depth-channel MRFT system identification and PID tuning workflow
%
% Purpose:
%   1. Convert experimental MRFT oscillation data into complex FRF samples.
%   2. Fit the selected low-order depth model in the Nyquist domain.
%   3. Build the identified transfer function and state-space realization.
%   4. Tune and compare several PID controllers:
%        - GA/ISE-based PID using the identified model
%        - MRFT GM-beta tuner
%        - Ziegler-Nichols from MRFT
%        - Tyreus-Luyben from MRFT
%
% Notes:
%   - Helper functions are assumed to be available on the MATLAB path.
%   - The numerical/code lines are intentionally preserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;clc;close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Experimental MRFT data: depth channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each row corresponds to one MRFT run:
%   column 1: beta value
%   column 2: oscillation frequency [Hz]
%   column 3: oscillation amplitude

data_depth = [ ...
 0.5   0.033908046  0.396075;
 0.4   0.043478261  0.29585;
 0.3   0.051315789  0.225525;
 0.2   0.058823529  0.18767775;
 0.1   0.064583333  0.158525;
 0     0.076923077  0.122225;
-0.1   0.083333333  0.106875;
-0.2   0.087121212  0.097925;
-0.3   0.1          0.074925;
-0.4   0.1125       0.070325;
-0.5   0.118055556  0.06085];


% Separate the experimental data into vectors used by the identification
% routine. Frequency is converted from Hz to rad/s for frequency-domain use.

beta_exp_vec    = data_depth(:,1);
freqHz_exp_vec  = data_depth(:,2);
amp_exp_vec     = data_depth(:,3);
freqOm_exp_vec  = 2*pi*freqHz_exp_vec;   % rad/s


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Convert MRFT oscillation data into complex FRF samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The describing-function relation maps each MRFT run into one complex
% Nyquist point. These points are the frequency-domain identification data.

Nb         = numel(beta_exp_vec);
W_exp_pts  = NaN(1, Nb);
mrft_h     = 5;
h_est      = mrft_h;

for i_B = 1:Nb
    beta_used       = beta_exp_vec(i_B);
    Amp             = amp_exp_vec(i_B);
    W_exp_pts(i_B)  = -(pi*Amp/(4*h_est)) * ...
                      ( sqrt(1-beta_used^2) + 1j*beta_used );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Select MRFT points used for fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By default, all available MRFT samples are used. If some runs are outliers
% or should be excluded, modify include_idx only.

include_idx = 1:Nb;      % use all samples by default

w_inc       = freqOm_exp_vec(include_idx);
Wd_inc      = W_exp_pts(include_idx).';

w  = w_inc(:);           % column vectors
Wd = Wd_inc(:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Select and initialize the system-identification model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is configured for the case_a depth-channel structure.

model_type = 'case_a';
fprintf('Using model: %s\n', upper(model_type));

[Wmod, th0] = create_caseA_sysid_model(w);   % local helper below

mag_floor = 1e-6;                            % avoid division by 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Identification weighting strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wts01 controls the relative influence of each MRFT frequency point in the
% fitting cost. Unity weights treat all selected points equally.

wts01 = ones(numel(w),1);  % default unity weights

% Example alternatives (uncomment ONE if needed and comment the others):
% wts01 = linspace(0.2,1,numel(beta_exp_vec(include_idx))).';
% wts01 = linspace(1,1,numel(beta_exp_vec(include_idx))).';
% wts01 = gaussian_beta_weights(beta_exp_vec(include_idx), 1, 0.2, true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Candidate cost functions for complex-domain fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Several fitting criteria are kept here for quick switching. MATLAB uses
% the last active assignment to cost.

cost = @(th) mean( abs( (Wmod(th,w)-Wd)./(abs(Wmod(th,w))+mag_floor) ).^2 );

cost = @(th) mean(  abs( (Wmod(th,w)-Wd)./(abs(Wd)+mag_floor) ).^2 );

cost = @(th) mean( ( log(abs(Wmod(th,w)./Wd))).^2 + ...
                    (angle(Wmod(th,w)) - angle(Wd)).^2 );

cost = @(th) mean( wts01 .* abs( (Wmod(th,w)-Wd) ./ ...
                    (abs(Wd)+mag_floor) ).^2 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. System identification using GA followed by local polishing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The initial parameter vector is scaled to define broad search bounds.
% GA performs global exploration; fminsearch then locally polishes the fit.

g0    = th0(:).';
scale = max(1, abs(g0));
lb    = -100*scale;
ub    =  1000*scale;
ub(end) = 2;0.3;   % delay upper bound

[th_hat, fval] = ga_sysid_search(cost, g0, ...
    'lb', lb, 'ub', ub, ...
    'PopulationSize', 5e4, 'MaxGenerations', 100, ...
    'UseParallel', false, 'Polish', false);

opts   = optimset('Display','iter','MaxFunEvals',5e4,'MaxIter',5e4);
% opts = optimset('Display','iter','MaxFunEvals',5e4,'MaxIter',5e4, 'TolX', 1e-8, 'TolFun', 1e-8);

th_hat = fminsearch(cost, th_hat, opts);     % local polish

Wfit = Wmod(th_hat,w);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8. Nyquist comparison: MRFT-derived data vs fitted model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This plot gives the primary visual check of the frequency-domain fit.

figure;
plot(real(Wd),  imag(Wd), 'o', ...
     real(Wfit),imag(Wfit),'x', 'LineWidth',1.2);
grid on; axis equal;
legend('MRFT FRF','Fitted model');
xlabel('Re'); ylabel('Im');
title('Nyquist: data vs fit');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 9. Display identified model information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print_sysid_model provides a formatted summary if available on the path.

print_sysid_model(model_type, th_hat, w);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 10. Map the fitted parameter vector to physical model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The fitted vector is converted into the positive parameter values used in
% the depth-channel transfer-function structure.

k1     = abs(th_hat(1));
Th     = abs(th_hat(2));
Tth    = abs(th_hat(3));
zeta   = abs(th_hat(4));
tau_id = max(0, th_hat(5));

fprintf('\nIdentified parameter vector th_hat:\n');
disp(th_hat(:).');

fprintf('case_a mapping:\n');
fprintf('  k1   = %.6g\n', k1);
fprintf('  Th   = %.6g\n', Th);
fprintf('  Tth  = %.6g\n', Tth);
fprintf('  zeta = %.6g\n', zeta);
fprintf('  tau  = %.6g\n\n', tau_id);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 11. Build identified depth transfer function and state-space model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The delay-free part of the identified model is converted to state space.
% The delay tau_id is retained separately for CLTools simulations.

s    = tf('s');

den1 = [1 0];                       % s
den2 = [Th 1];                      % (Th s + 1)
den3 = [Tth^2, 2*zeta*Tth, 1];      % (Tth^2 s^2 + 2 ζ Tth s + 1)

G_id = tf(k1, conv(conv(den1, den2), den3));

G_trunc        = truncate_tf_coeffs(G_id, 1e-2);
[Ai,Bi,Ci,Di]  = ssdata(ss(G_trunc));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 12. Closed-loop simulation settings for PID tuning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These options define the common simulation environment used by all PID
% comparison methods.

dt   = 1e-3;       % time step (should be < tau_id)
Tsim = 50;         % simulation horizon [s]
u_f  = 1;          % unconstrained factor (1 = constrained)

opts = struct( ...
    'PadeOrder', 4, ...
    'ClampU', true, 'UMax', 200*u_f, 'UMin', -200*u_f, ...
    'DFilter', true, 'DTau', 0.05, ...
    'AntiWindup', true, 'Imin', -100*u_f, 'Imax', 100*u_f, ...
    'Beta', 1 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 13. PID tuning using GA/ISE on the identified model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The PID gains are optimized directly on the identified state-space model
% using the CLTools closed-loop simulation cost.

pid_cost = @(g) CLTools.sim_cost(Ai,Bi,Ci,Di, tau_id, ...
                abs(g(1)), abs(g(2)), abs(g(3)), dt, Tsim, opts);

% starting guess and bounds
g0 = [33.6564143666146, 1.05002086367551e-10, 81.6719991685690];
lb = [0 0 0];
ub = [1000 1000 500];

[gopt, fval_pid] = ga_pid_search(pid_cost, g0, ...
    'lb', lb, 'ub', ub, ...
    'PopulationSize', 100, 'MaxGenerations', 50, ...
    'UseParallel', false, 'Polish', false);

Kp_ISE = gopt(1);
Ki_ISE = abs(gopt(2));
Kd_ISE = gopt(3);

fprintf('GA ISE PID: J=%.6g | Kp=%g Ki=%g Kd=%g\n', ...
        fval_pid, Kp_ISE, Ki_ISE, Kd_ISE);

Jcheck = CLTools.sim_cost(Ai,Bi,Ci,Di, tau_id, ...
                          Kp_ISE, Ki_ISE, Kd_ISE, dt, Tsim, opts);

fprintf('Direct cost (ISE, GA gains): %.6g\n', Jcheck);

[ti, yi, ~] = CLTools.sim_cl(Ai,Bi,Ci,Di, tau_id, ...
                             Kp_ISE, Ki_ISE, Kd_ISE, dt, Tsim, opts);

figure;
plot(ti, yi, 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('y');
title(sprintf('PID-B tuned by ISE (GA): Kp=%.3g, Ki=%.3g, Kd=%.3g', ...
              Kp_ISE, Ki_ISE, Kd_ISE));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 14. PID tuning using the MRFT GM-beta method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This tuner uses the experimental MRFT data directly, without using the
% fitted transfer function as the tuning basis.

beta_set = -0.2; GM_set = 2;

[Kp_GMB, Ki_GMB, Kd_GMB, ~, ~, ~] = PIDgains_GM_Beta( ...
    mrft_h, beta_exp_vec, amp_exp_vec, freqHz_exp_vec, ...
    'beta', beta_set, 'GM', GM_set);

fprintf('GM–beta PID: Kp=%.4g, Ki=%.4g, Kd=%.4g\n', ...
        Kp_GMB, Ki_GMB, Kd_GMB);

[tg, yg, ~] = CLTools.sim_cl(Ai,Bi,Ci,Di, tau_id, ...
                             Kp_GMB, Ki_GMB, Kd_GMB, dt, Tsim, opts);

figure;
plot(tg, yg, 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('y');
title(sprintf('GM–beta PID step: Kp=%.3g, Ki=%.3g, Kd=%.3g', ...
              Kp_GMB, Ki_GMB, Kd_GMB));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 15. PID tuning using Ziegler-Nichols from MRFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The ultimate gain and period are inferred from the MRFT dataset and then
% converted to PID gains using the ZN rule.

[zn, Kp_ZN, Ki_ZN, Kd_ZN] = zn_from_mrft( ...
    beta_exp_vec, amp_exp_vec, freqHz_exp_vec, mrft_h, 'PID');

fprintf('ZN from MRFT (beta≈%g): Ku=%.4g, Tu=%.4g s\n', ...
        zn.beta_used, zn.Ku, zn.Tu);

fprintf('ZN PID: Kp=%.4g, Ki=%.4g, Kd=%.4g\n', ...
        Kp_ZN, Ki_ZN, Kd_ZN);

[tz, yz, ~] = CLTools.sim_cl(Ai,Bi,Ci,Di, tau_id, ...
                             Kp_ZN, Ki_ZN, Kd_ZN, dt, Tsim, opts);

figure;
plot(tz, yz, 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('y');
title(sprintf('ZN-from-MRFT PID step: Kp=%.3g, Ki=%.3g, Kd=%.3g', ...
              Kp_ZN, Ki_ZN, Kd_ZN));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 16. PID tuning using Tyreus-Luyben from MRFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The same MRFT-derived ultimate quantities are used with the TL tuning rule,
% which typically gives more conservative PID gains than ZN.

[tl, Kp_TL, Ki_TL, Kd_TL] = tl_from_mrft( ...
    beta_exp_vec, amp_exp_vec, freqHz_exp_vec, mrft_h, 'PID');

fprintf('TL from MRFT (beta≈%g): Ku=%.4g, Tu=%.4g s\n', ...
        tl.beta_used, tl.Ku, tl.Tu);

fprintf('TL PID: Kp=%.4g, Ki=%.4g, Kd=%.4g\n', ...
        Kp_TL, Ki_TL, Kd_TL);

[tt, yt, ~] = CLTools.sim_cl(Ai,Bi,Ci,Di, tau_id, ...
                             Kp_TL, Ki_TL, Kd_TL, dt, Tsim, opts);

figure;
plot(tt, yt, 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('y');
title(sprintf('TL-from-MRFT PID step: Kp=%.3g, Ki=%.3g, Kd=%.3g', ...
              Kp_TL, Ki_TL, Kd_TL));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End of script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%