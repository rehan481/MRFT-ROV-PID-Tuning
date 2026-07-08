function G_trunc = truncate_tf_coeffs(G_id, rel_tol)
% TRUNCATE_TF_COEFFS Truncates relatively small coefficients in the transfer function and normalizes to monic denominator.
%   G_trunc = truncate_tf_coeffs(G_id, rel_tol) takes a transfer function G_id
%   and a relative tolerance rel_tol (default 1e-8), sets coefficients smaller
%   than rel_tol times the maximum absolute coefficient (in num and den separately,
%   after scaling each to max abs=1) to zero, removes leading zeros, normalizes
%   the denominator to have leading coefficient 1 (by dividing num and den by the
%   leading den coefficient), and returns the truncated and normalized transfer function G_trunc.

if nargin < 2
    rel_tol = 1e-8;
end

[num, den] = tfdata(G_id, 'v');

% Handle zero numerator or denominator cases
scale_num = max(abs(num));
if scale_num == 0
    G_trunc = tf(0);
    return;
end

scale_den = max(abs(den));
if scale_den == 0
    error('Denominator is zero');
end

% Scale to max abs=1
num_scaled = num / scale_num;
den_scaled = den / scale_den;

% Set relatively small coefficients to zero
num_scaled(abs(num_scaled) < rel_tol) = 0;
den_scaled(abs(den_scaled) < rel_tol) = 0;

% Remove leading zeros
while ~isempty(num_scaled) && num_scaled(1) == 0
    num_scaled(1) = [];
end
if isempty(num_scaled)
    num_scaled = 0;
end

while ~isempty(den_scaled) && den_scaled(1) == 0
    den_scaled(1) = [];
end
if isempty(den_scaled) || den_scaled(1) == 0
    error('Denominator became zero after truncation');
end

% Reconstruct with original scales
num_trunc = num_scaled * scale_num;
den_trunc = den_scaled * scale_den;

% Normalize to monic denominator (divide num and den by leading den coeff)
if ~isempty(den_trunc) && den_trunc(1) ~= 0
    factor = den_trunc(1);
    num_trunc = num_trunc / factor;
    den_trunc = den_trunc / factor;
end

% Create the truncated and normalized transfer function
G_trunc = tf(num_trunc, den_trunc);

end