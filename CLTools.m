classdef CLTools
    % Minimal closed-loop helpers (PID Type-B, Euler, ISE)

    methods (Static)


        function J = sim_cost(A,B,C,D,tau,Kp,Ki,Kd,dt,T, varargin)
            [t,y] = CLTools.sim_cl(A,B,C,D,tau, abs(Kp),abs(Ki),abs(Kd), dt,T, varargin{:});


            %             [~,y] = CLTools.sim_cl(A,B,C,D,tau,Kp,Ki,Kd,dt,T);
            %%% OR

            %                     [~,y] = CLTools.sim_cl(A,B,C,D,tau, abs(Kp),abs(Ki),abs(Kd), dt,T, ...
            %                         'ClampU',true,'UMax',200,'UMin',-200, ... % ',1650,'UMin',1350, ..
            %                         'DFilter',true,'DTau',0.05, ...
            %                         'AntiWindup',true,'Imin',-100,'Imax',100, ...
            %                         'Beta',1);

            e = 1 - y;
%             J = sum(e.^2)*dt;   % ISE
%             J = trapz(t, t .* abs(e));   % ITAE
            J = 10*sum(e.^2)*dt+trapz(t, t .* abs(e));   % ITAE
        end


      
        function [t,y,u] = sim_cl(A,B,C,D,tau,Kp,Ki,Kd,dt,T, varargin)
            % Using Pade approximation for the Delay. FOr orignal Delay
            % implementation  use the commented code below (commnt ths Pade
            % verion and Uncmnt below version)

            opts = CLTools.parse_opts(varargin{:});


            % -------- init ----------
            N  = round(T/dt)+1;  t = (0:N-1)'*dt;
            nx = size(A,1);      x = zeros(nx,1);
            yk = 0;  u = 0;      % current y and u (scalars)
            y_prev = 0;          % for dy
            have_prev_y = false;

            % Augment system with Pade approximation if specified
            if opts.PadeOrder > 0
                order = opts.PadeOrder;
                [num, den] = pade(tau, order);
                sys_delay = tf(num, den);
                sys_plant = ss(A, B, C, D);
                sys_total = series(sys_delay, sys_plant);
                A = sys_total.A;
                B = sys_total.B;
                C = sys_total.C;
                D = sys_total.D;
                nx = size(A,1);
                x = zeros(nx,1);
                dN = 0;  % no discrete delay buffer needed
            else
                dN = max(0, round(tau/dt));  % discrete pure delay
            end

            ub = zeros(dN+1,1);

            yvec = zeros(N,1);
            uvec = zeros(N,1);

            % derivative filter state (if enabled)
            ydot_f_prev = 0;
            ydot_f_initialized = false;

            % integrator state
            ui = 0;

            % -------- loop ----------
            for k = 1:N
                r = 1;                % unit-step reference
                e = r - yk;           % error

                % --- Proportional (Type-B) ---
                p_term = Kp * (opts.Beta * r - yk);

                % --- Integral ---
                ui = ui + e*dt;
                if opts.AntiWindup
                    % clamp integrator state
                    if ui > opts.Imax, ui = opts.Imax; end
                    if ui < opts.Imin, ui = opts.Imin; end
                end
                i_term = Ki * ui;

                % --- Derivative on measurement (with optional 1st-order filter) ---
                if k == 1
                    dy_meas = 0;   % first tick: no derivative yet
                else
                    dy_meas = (yk - y_prev)/dt;
                end

                if opts.DFilter && opts.DTau > 0
                    alpha = dt / (opts.DTau + dt);         % same as in your C++
                    if ~ydot_f_initialized
                        ydot_f_prev = dy_meas;
                        ydot_f_initialized = true;
                    end
                    ydot_used = ydot_f_prev + alpha*(dy_meas - ydot_f_prev);
                    ydot_f_prev = ydot_used;
                else
                    ydot_used = dy_meas;
                end

                d_term = -Kd * ydot_used;  % D on measurement ⇒ minus sign

                % --- raw control and optional clamping ---
                u_raw = p_term + i_term + d_term;
                if opts.ClampU
                    u = min(max(u_raw, opts.UMin), opts.UMax);
                else
                    u = u_raw;
                end

                % --- handle pure delay on control (if applicable) ---
                if dN > 0
                    ub = [u; ub(1:end-1)];
                    ud = ub(end);
                else
                    ud = u;
                end

                % --- plant Euler step ---
                x = x + dt*(A*x + B*ud);
                y_prev = yk;
                yk = C*x + D*ud;

                % store
                yvec(k) = yk;
                uvec(k) = u;
            end

            y = yvec;
            u = uvec;
        end

        % -------- option parsing helper --------
        function opts = parse_opts(varargin)
            % defaults (match legacy behavior)
            opts.ClampU     = false;
            opts.UMin       = -Inf;
            opts.UMax       = +Inf;
            opts.DFilter    = false;
            opts.DTau       = 0.05;
            opts.AntiWindup = false;
            opts.Imin       = -100;
            opts.Imax       = +100;
            opts.Beta       = 1.0;
            opts.PadeOrder  = 0;
            opts.Bd = [];   % optional custom disturbance mapping (nx x 1). If empty, use B.


            % NEW: constant plant disturbance & optional static trim
            opts.BuoyancyN  = 0;   % constant buoyancy/disturbance (units = model input)
            opts.UTrim      = 0;   % constant controller output bias (input units)

            % ---- struct form ----
            if nargin == 1 && isstruct(varargin{1})
                s = varargin{1};
                f = fieldnames(s);
                for i = 1:numel(f)
                    if ~isfield(opts, f{i})
                        error('Unknown option "%s".', f{i});
                    end
                    opts.(f{i}) = s.(f{i});
                end
            else
                % ---- name-value pairs ----
                if mod(nargin,2) ~= 0
                    error('Optional arguments must be name-value pairs or a struct.');
                end
                for i = 1:2:nargin
                    name = varargin{i};
                    val  = varargin{i+1};
                    if ~ischar(name) && ~isstring(name)
                        error('Option name must be char or string.');
                    end
                    name = char(name);
                    if ~isfield(opts, name)
                        error('Unknown option "%s".', name);
                    end
                    opts.(name) = val;
                end
            end

            % ---- sanity ----
            if ~(opts.UMin <= opts.UMax), error('UMin must be <= UMax.'); end
            if ~(opts.Imin <= opts.Imax), error('Imin must be <= Imax.'); end
            if opts.DTau < 0, error('DTau must be >= 0.'); end
            if opts.PadeOrder < 0 || abs(round(opts.PadeOrder)-opts.PadeOrder) > eps
                error('PadeOrder must be a nonnegative integer.');
            end
            if ~isfinite(opts.BuoyancyN), error('BuoyancyN must be finite.'); end
            if ~isfinite(opts.UTrim),     error('UTrim must be finite.'); end
        end







    end

end
