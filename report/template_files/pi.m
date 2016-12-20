% Generate random points inside a centered square of area 4
trials = 1e6;
x = 2*rand(trials,1) - 1;
y = 2*rand(trials,1) - 1;

% Distance of points to (0,0)
r = sqrt(x.^2 + y.^2);

% Fraction of points inside unit circle
frac_inside = sum(r < 1)/trials;

fprintf('Pi is approximately %.6f.\n', 4*frac_inside);
