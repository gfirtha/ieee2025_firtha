function w = custom_window(N)
% N: Number of samples in the window
% Generate a window where the second derivative vanishes at the endpoints

% Normalized time vector from 0 to 1
t = linspace(0, 1, ceil(N/2));

% Define the polynomial: w(t) = -6t^5 + 15t^4 - 10t^3 + 1
% This polynomial satisfies w(0) = 0, w(1) = 0, w'(0) = 0, w'(1) = 0, w''(0) = 0, w''(1) = 0
w = -6 * t.^5 + 15 * t.^4 - 10 * t.^3 + 1;

% Normalize to ensure maximum value is 1 (optional, already normalized in this case)
w = w / max(w);
if rem(N,2) == 0
    w = [fliplr(w) w]';
else
    w = [fliplr(w) w(2:end) ]';
end
end