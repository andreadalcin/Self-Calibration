function focal_lengths = getInitialEstimateOfFocalLength(F, width, height, opening_angles, draw_plot)

tau = 1;     % Aspect ratio = 1

% Principal points
u0 = width / 2;
v0 = height / 2;

focal_lengths = zeros(size(opening_angles,2), 1);

for i = 1:size(opening_angles,2)
    f0 = max(width, height) / 2 / tan(deg2rad(opening_angles(i) / 2));

    % Undo intrinsic parameters
    G = [tau, 0, 0; 0, 1, 0; u0, v0, 1] * F * [tau, 0, u0; 0, 1, v0; 0, 0, 1];
    G_1 = diag([f0,f0,1]) * G * diag([f0,f0,1]);
    G_1 = G_1 / norm(G_1, 'fro');

    % SVD
    [U,S,V] = svd(G_1);
    u13 = U(1,3);
    u23 = U(2,3);
    v13 = V(1,3);
    v23 = V(2,3);
    a = S(1,1);
    b = S(2,2);

    % Quadratic equation
    c1 = a^2 * (1 - u13^2) * (1 - v13^2) - b^2 * (1 - u23^2) * (1 - v23^2);
    c2 = a^2 * (u13^2 + v13^2 - 2 * u13^2 * v13^2) - b^2 * (u23^2 + v23^2 - 2 * u23^2);
    c3 = a^2 * u13^2 * v13^2 - b^2 * u23^2 * v23^2;

    roots = zeros(2,1);
    d = sqrt(c2^2 - 4*c1*c3);
    roots(1) = ( -c2 + d ) / (2*c1);
    roots(2) = ( -c2 - d ) / (2*c1);

    % De-normalize
    roots = sqrt(roots) * f0;
    [~,I] = min(abs(roots - f0)); % Add linear equations check?
    focal_lengths(i) = roots(I);
end

if draw_plot
    figure;
    plot(real(focal_lengths), '.');
end

focal_lengths = nonzeros(real(focal_lengths));
focal_lengths = focal_lengths(~isnan(focal_lengths));
end