from sys import argv

eps = float(argv[1])

real_len = 1.0
N = 200
M = 55

dx = real_len / (N + eps)
solid_len = (M + 1 - eps) * dx
total_len = real_len + solid_len

print(
    "geometry.prob_hi=%20.16f\t%20.16f\t%20.16f   #vfrac %5.3f"
    % (total_len, total_len / 16.0, total_len / 16.0, eps)
)
