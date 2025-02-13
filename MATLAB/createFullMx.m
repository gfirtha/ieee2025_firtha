function [xOut,yOut, B] = createFullMx(x,y,A)

% Suppose A is your original N x N matrix

% Create the top-left quadrant (mirror in both x and y, but skip the first row and column to avoid duplication)
A_corner = flipud(fliplr(A(2:end,2:end)));

% Create the top-right quadrant (mirror in y only)
A_top = flipud(A(2:end,:));

% Create the bottom-left quadrant (mirror in x only)
A_left = fliplr(A(:,2:end));

% The bottom-right quadrant is just the original A.
% Now, stitch them together:
B = [ A_corner, A_top;
      A_left,   A ];

xOut = [flip(-x(2:end));x];
yOut = [flip(-y(2:end));y];

end

