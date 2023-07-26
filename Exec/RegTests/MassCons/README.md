# MassCons

This is a simple contrived test case to demonstrate mass conservation. A cubical
domain with various types of boundary (NoSlipWall, SlipWall, and Symmetry) is
initialized with a uniform velocity field with components in all directions. It
is run with all of the major numerical schemes in PeleC (Godunov/PLM,
Godunov/PPM, MOL/PLM, MOL/PCM). The case is nonreacting with GammaLaw EOS, but
diffusion is turned on. This case can also test implementation of walls through
bcnormal when BCs are set to "Hard". In that case, wall_type = 0 or 1 are
NoSlip and Slip walls, respectively.
