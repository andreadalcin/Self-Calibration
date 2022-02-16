function f0 = getInitialEstimateOfFocalLength(F, width, height)
f0 = 2 * (width + height);
f0 = getInitialEstimateFromF0(F, f0, width, height);
% focal_lengths = nonzeros(real(focal_lengths));
% focal_lengths = focal_lengths(~isnan(focal_lengths));
end