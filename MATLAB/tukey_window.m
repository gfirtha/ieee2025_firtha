function w = tukey_window(L, alpha)
% TUKEY_WINDOW Compute Tukey window of length L with parameter alpha
%   w = tukey_window(L, alpha) returns a vector w of length L containing
%   the Tukey window samples. alpha is the taper parameter (0 <= alpha <= 1).
%   Does not use MATLAB's built-in tukeywin function.

% Input validation
if L < 1 || floor(L) ~= L
    error('L must be a positive integer');
end
if alpha < 0 || alpha > 1
    error('alpha must be between 0 and 1');
end

% Initialize output vector
w = zeros(L, 1);

% Handle alpha = 0 (rectangular window)
if alpha == 0
    w = ones(L, 1);
    return;
end

% Time indices (0 to L-1, corresponding to t from 0 to L in continuous time)
n = (0:L-1)';

% Normalized time: map n to t in [0, L]
t = n * (L / (L-1)); % Adjust to span [0, L] exactly

% Compute boundaries of the tapered regions
taper_length = alpha * L / 2;

% Left taper: 0 <= t <= alpha*L/2
left_mask = t <= taper_length;
w(left_mask) = 0.5 * (1 + cos((2*pi/(alpha*L)) * (t(left_mask) - alpha*L/2)));

% Middle flat region: alpha*L/2 < t < L - alpha*L/2
middle_mask = (t > taper_length) & (t < L - taper_length);
w(middle_mask) = 1;

% Right taper: L - alpha*L/2 <= t <= L
right_mask = t >= L - taper_length;
w(right_mask) = 0.5 * (1 + cos((2*pi/(alpha*L)) * (t(right_mask) - L + alpha*L/2)));

end